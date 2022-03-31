
#!/usr/bin/perl
#author:wangxin
### function: To determine the gene existence status between mutation samples and  non-mutation samples,
### those samples are treated with immunotherapy
### 
### here is to check whether gene mutations when treated with immnunotherapy can have a better survival

### 

use strict;
use warnings;
use Statistics::Test::WilcoxonRankSum;
use List::Util qw(sum);
use List::Util 'shuffle';
use Statistics::R;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","g:s","o:s","t:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i}|| !defined $opts{g} ||!defined $opts{t}||!defined $opts{o}) {
       	die "************************************************************************
       	Usage: $0.pl
	       		-i: mutation information of each gene from each samples 
			-t: specfic gene list for BRCA genetic instability
			-g: clinical data
			-o: output of pvalue for survival analysis
************************************************************************\n";
}



my $index=$opts{i};
my $clinic=$opts{g};
my $Gene_ID=$opts{t};


######index file to calculate whether they have gene mutation for the specific samples 

my %hash;my $n; my %mutant;
open INDEX,"$index" or die $!;
while (<INDEX>){
	chomp;
	s/\r//;
    my ($gene,$Variant_ty,$sample)=(split/\t/,$_)[0,9,16];
	next if ($gene eq "Hugo_Symbol");
	
	my $sp= join "-",(split/\-/,$sample)[0,1];
	$hash{$sp}++;
	if ($Variant_ty eq "Nonsense_Mutation" || ($Variant_ty eq "Missense_Mutation" ) || $Variant_ty eq "In_Frame_Del" ||$Variant_ty eq "Frame_Shift_Del" || $Variant_ty eq "Frame_Shift_Ins" || $Variant_ty eq "In_Frame_Ins" ||$Variant_ty eq "Nonstop_Mutation" || $Variant_ty eq "Translation_Start_Site" || $Variant_ty eq  "Splice_Site"){
		$mutant{$gene}->{$sp}++;
	}
}

close INDEX;

# my @sample;
# foreach my $i ( sort {$hash{$b} <=> $hash{$a}} keys %hash ){
# 	push @sample,$i;
#
# }

###### index for the clinical survival information
my %time_s; my %status_s; my %ml; 
open CL, "$clinic" or die $!;
while (<CL>){
	chomp;
	next if ($_=~/^#/);
	next if ($_=~/^OTHER/);
	my ($sam,$mutationload,$ttime,$sstatus)=(split/\t/,$_)[0,1,3,4];
	#print "$ttime\n";
	$time_s{$sam}=$ttime;
	$ml{$sam}=$mutationload;
	#next if ($sstatus eq "Not")
	next if ($sstatus =~/Available/g || $ttime =~ /Available/g);
	$status_s{$sam}= ($sstatus eq "LIVING")?0:1;
}

close CL;

my $out=$opts{o};
open STAT,">$out" or die $!;
print STAT "Gene\tPvalue_mutationload\tPvalue_survival\tMutationload_up\tMutationload_down\tSurvival_up\tSurvival_down\n";

#### reading the potential gene markers for the patients (Here the gene markers are the genes that are identified from )
open GL,"$Gene_ID" or die $!;

while (<GL>){
	chomp;
	s/\r//;
	my $ge=$_;
	
	
	#### input of samples then divided them into two categories: gene mutants and gene nonmutants
	open OUT,">tmp.txt" or die $!;
	
	print OUT "time\tstatus\tx\tMutationload\n";
	
	my @time_up; my @time_down; my @mut_up; my @mut_down;
	foreach my $i (keys %hash){
		if (exists $mutant{$ge}->{$i}){
			print OUT "$time_s{$i}\t$status_s{$i}\tMutant\t$ml{$i}\n";
			push @time_up,$time_s{$i};
			push @mut_up,$ml{$i};
		}else{
			
			print OUT"$time_s{$i}\t$status_s{$i}\tNonMutant\t$ml{$i}\n";
			push @time_down,$time_s{$i};
			push @mut_down,$ml{$i};
		}
	}
	
	# patients number
	
	my $size_up=scalar @mut_up;
	my $size_down=scalar @mut_down;
	
	
	if (scalar  @mut_up <=1 ){
		print STAT "$ge\tNA\tNA\tNA\tNA\tNA\tNA\t$size_up\t$size_down\n";
		next;
	}
	
	
	my $MLaverage_up=mean(@mut_up);
	my $MLaverage_down=mean(@mut_down);
	
	my $wilcoxon_test= Statistics::Test::WilcoxonRankSum->new();
	
	$wilcoxon_test->load_data(\@mut_up, \@mut_down);
	my $pf= sprintf '%f',$wilcoxon_test->probability();
	
	
	my $average_up=mean(@time_up);
	my $average_down=mean(@time_down);
	

	
	close OUT;
	
	system ("cp tmp.txt $ge.txt");
	my $R=Statistics::R->new();
	$R->startR ;
	$R -> run('library(survival);');
	
	# $R->set('x',\@numbers)

	$R -> run('table <- read.table("tmp.txt", header=T);');
	$R -> run('attach(table);');
	$R -> run('fit<-survdiff(Surv(time, status)~x, data=table);');
	
	$R -> run('pvalue<- 1 - pchisq(fit$chisq, length(fit$n) - 1)');
	##
	#$R -> run('pvalue<-ifelse ( is.na(fit),next,(round(1 - pchisq(fit$chisq, length(fit$n) - 1),3)))[[1]]');
	#$R -> send('print(meang1);');

	my $pvalue= $R->get('pvalue');
	
	print STAT "$ge\t$pf\t$pvalue\t$MLaverage_up\t$MLaverage_down\t$average_up\t$average_down\t$size_up\t$size_down\n";
	$R->stopR();
	system ("rm tmp.txt");
	
}


close STAT;

#
sub mean {
	return sum(@_)/@_;
}
