#!/usr/bin/perl


# author:wangxin
# contact: xin.wang@childrens.harvard.edu
# PI: kaifu Chen
# contact: Kaifu.chen@childrens.harvard.edu

#######################################################################################################################################################################
#### Function: To extract TCGA samples with radiations 
####			 
#######################################################################################################################################################################

use strict;
use warnings;
use Statistics::Test::WilcoxonRankSum;
use Statistics::TTest;
use List::Util qw(sum min max shuffle);
use Statistics::Descriptive;
use Statistics::Basic qw(:all);
use List::Uniq ':all';

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","g:s","s:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i}|| !defined $opts{g} || !defined $opts{s} ) {
       	die "************************************************************************
       	Usage: $0.pl
	       		-i: query files of nonsilent mutations statistics (Mutation_information.txt)
			-g: index files of sample mutation load (mc3)
			-s: cancer type
************************************************************************\n";
}


my $input=$opts{i};
my %hash;
my $index=$opts{g};
my $cancer=$opts{s};
#my $potential=$opts{p};

#### build the mutation load index of each patient under required cancer type.
my %mutation;
open I,"$input" or die $!;
while (<I>){
	chomp;
	s/\r//;
	my ($sample,$type,$all)=(split/\t/,$_)[0,1,2];
	next unless ($type eq $cancer);	
	### here we ignore the patient that have no mutations at all
	next if ($all<=0);
	$mutation{$sample}=$all;
	#print "$type\t$mutation{$type}->{$sample}\n";
}
close I;


#### build the deleterious index of each gene in each patient.
my %nonmut; my %type2; my %geneindex; my %bad; my %nondel; my %Sampbad;
open INDEX,"$index" or die $!;
while (<INDEX>){
	chomp;

	s/\r//g;
	my ($gene,$Variant_ty,$ref,$alt,$sample,$deleterious,$damage)=(split/\t/,$_)[0,8,11,12,15,71,72];
	next if ($gene eq "Hugo_Symbol");
	my $Var_change=join "/",($ref,$alt);
	my $sp= join "-",(split/\-/,$sample)[0,1,2,3];
	$sp=~s/[A-Z]$//;
	next unless (exists $mutation{$sp});
	$nonmut{$gene}->{$sp}++;
	$hash{$gene}->{$sp}->{$Variant_ty}++;
	
	
	#$type2{deleterious}++ if ($deleterious =~/deleterious/ || $damage =~ /probably_damaging/ || $Variant_ty eq "Nonsense_Mutation" || $Variant_ty  eq "Frame_Shift_Del" || $Variant_ty eq "Splice_Site" || $Variant_ty eq "Frame_Shift_Ins" || $Variant_ty eq "Nonstop_Mutation");

	
	#### index of deleterious mutation samples
	$type2{deleterious}++ if ($deleterious =~/deleterious/ || $damage =~ /probably_damaging/ || $Variant_ty eq "Nonsense_Mutation" || $Variant_ty  eq "Frame_Shift_Del" || $Variant_ty eq "Splice_Site" || $Variant_ty eq "Frame_Shift_Ins" || $Variant_ty eq "Nonstop_Mutation" || $Variant_ty  eq "In_Frame_Del" ||$Variant_ty  eq "In_Frame_Ins" );
	$bad{$gene}->{$sp}++ if ($deleterious =~/deleterious/ || $damage =~ /probably_damaging/ || $Variant_ty eq "Nonsense_Mutation" || $Variant_ty  eq "Frame_Shift_Del" || $Variant_ty eq "Splice_Site" || $Variant_ty eq "Frame_Shift_Ins" || $Variant_ty eq "Nonstop_Mutation"|| $Variant_ty  eq "In_Frame_Del" ||$Variant_ty  eq "In_Frame_Ins");
	
	$Sampbad{$sp}->{$gene}++ if ($deleterious =~/deleterious/ || $damage =~ /probably_damaging/ || $Variant_ty eq "Nonsense_Mutation" || $Variant_ty  eq "Frame_Shift_Del" || $Variant_ty eq "Splice_Site" || $Variant_ty eq "Frame_Shift_Ins" || $Variant_ty eq "Nonstop_Mutation"|| $Variant_ty  eq "In_Frame_Del" ||$Variant_ty  eq "In_Frame_Ins");
	
	#$geneindex{$gene}++ if ($deleterious =~/deleterious/ || $damage =~ /probably_damaging/ || $Variant_ty eq "Nonsense_Mutation" || $Variant_ty eq "Frame_Shift_Del" || $Variant_ty eq "Splice_Site" || $Variant_ty eq "Frame_Shift_Ins" || $Variant_ty eq "Nonstop_Mutation"|| $Variant_ty  eq "In_Frame_Del" ||$Variant_ty  eq "In_Frame_Ins" );
	
	
	### index of non-deleteriou mutation samples
	#$nondel{$gene}->{$sp}++ if (($deleterious !~/deleterious/ && $damage !~ /probably_damaging/ && $Variant_ty eq "Missense_Mutation") || $Variant_ty eq "Silent" ||$Variant_ty eq "RNA" ||$Variant_ty eq "Intron" ||$Variant_ty eq "3'UTR" ||$Variant_ty eq "3'Flank" ||$Variant_ty eq "5'UTR" ||$Variant_ty eq "5'Flank");
	
	
	
}
close INDEX;




open OUT,">$cancer.mutationload.p1p2.txt" or die $!;


print OUT "Gene\tP1_wilcoin\tP2_wilcoin\tP1_ttest\tP2_ttest\tDel_mean\tCoDel_mean\tNonCoDel_mean\tDel_num\tCoDel_num\tNonCoDel_num\tDel_median\tCoDel_median\tNonCoDel_median\tDel_max\tCoDel_max\tNonCoDel_max\tDel_min\tCoDel_min\tNonCoDel_min\n";




#######################################################################################################################################################################
### Step 1. put all the deleterious mutants and mutation load into array;
###	We only consider those genes that have deleterious mutations.
#######################################################################################################################################################################

foreach my $n (keys %bad){
	#print  "$n\t";
	my @mut=(); my @non_mut=(); my @mut_name; my @non_mut_name; my $no=0; my $mu=0; my @silmut;  my $silmu=0;
	
	#### Step1. put all the mutation load of gene deleterious mutants into @mut;
	my %SM; my @SamMut=();
	## $i is sample $n is gene
	foreach my $i (keys %mutation){
		#print "$i\n";
		if(exists $bad{$n}->{$i}){
		push (@mut,$mutation{$i});
		push (@SamMut,$i);
		
		#### push index of sample into SM hash;
		$SM{$i}++;
		$mu++;
		}
	}
	
	##################################################################################################################################################################
	#### Step2. put mutation load of other mutants that did not show this gene deleterious but possess other co-deleterious mutations except this gene mutations. 
	### require: co-deleterious genes are those genes that share more than 60% of mutations in deleterious samples, 
	###          then other co-deleterious samples are those sample contain co-deleterious genes.
	##################################################################################################################################################################
	
	### Step2.1 find the other deleterious genes in this gene deleterious mutants. 
	
	my @OtherGenes=();
	my %numG;
	
	### $j means gene deleterious samples
	foreach my $j (@SamMut){
		
		### $q means all the deleterious genes in those gene deleterious samples
		foreach my $q (keys %{$Sampbad{$j}}){
			
			### ignore the deleterious genes and put all the co-mutated genes into hash $numG{$q}
			next if ($q eq $n);	
			push  @OtherGenes, $q;
			$numG{$q}++;
			### Test
			#print "$q\t";
		}
	}
	
	my @UniqOtherG=();
	
	
	#### We here further define the co-mutated genes that require muated in more than 60% of co-mutated samples, others are not co-mutated genes.
	foreach my $m (uniq(@OtherGenes)){
		if ($numG{$m} <=$mu && $numG{$m}>0.6*$mu ){
			push @UniqOtherG, $m;
		}
	}


	### we then use co-mutated genes and push other deleterious samples into array, ignoring the deleterious samples;
	
	my @Othedelmut=();
	my $Omut=0;
	my @SamOthers=();
	my %numS;
	
	### $g are the co-mutated genes
	foreach my $g (@UniqOtherG){
		
		### other samples is $h
		foreach my $h (keys %{$bad{$g}}){
			next if (exists $SM{$h});
			#print "$h\t";
			#push (@Othedelmut,$mutation{$h});
			push (@SamOthers,$h);
			
			$numS{$h}++;
			#$Omut++;
		}
	}
	
	
	### To improve the sensitivity of co-mutants, we require extral mutants cannot exceed 1.4* mut;
	
	### uniqOthersam are those samples with co-muations 
	my @uniqOtherSam=();
	foreach my $f (uniq(@SamOthers)){
		if ($numS{$f}< ($mu+0.4*$mu)){
			push @uniqOtherSam, $f;
		}
	}
	

	
	### calculate the number of co-muatation samples and put their mutation load into array and remember their id.
	my %target;
	foreach my $sam(@uniqOtherSam){
		next if (exists $SM{$sam});
		push (@Othedelmut,$mutation{$sam});
		$target{$sam}++;
		print "$mutation{$sam}\t";
		$Omut++;
	}
	
	print "\n";
	
	
	##################################################################################################################################################################
	#### Step3. put mutation load of other mutants that dont have gene and co-mutated genes deleterious mutation s 
	##################################################################################################################################################################
	
	### then push all the other samples without co-muation and without deleterious (Sampbad represent all the samples)
	my $Noncomut_num; my @Noncomut;
	foreach my $eachSam (keys %mutation){
		next if (exists $SM{$eachSam} || exists $target{$eachSam});
		push (@Noncomut,$mutation{$eachSam});	
		$Noncomut_num++;
	}
	


	#################################################################################################################################################################
	##### Step 4. Run the Wilcoxon test P1 to compare SamMut and Othedelmut, P2 to compare SamMut and Noncomut.
	#################################################################################################################################################################
	

	#### join the mutation load in the deleterious mutants and co-mutated mutants
	my $up_value= join "\t",@mut;
	#my $down_value2=join "\t",@non_mut;
	my $down_value=join "\t",@Othedelmut;
	

	#### build the wilcoxon test p1, here also added the t-test
	my $wilcoxon_test= Statistics::Test::WilcoxonRankSum->new();
    my $ttest=new Statistics::TTest;
	### eliminate those genes that only contain zero deleterious mutation (this might not possible) and zero co-mutated mutaitons 
	## 		if the gene is mutators, there will be more co-mutated mutations. So if the number is zero in co-muated genes, this would be the gene is not mutator.
	next unless (scalar @mut !=0 && scalar @Othedelmut !=0 && scalar @Noncomut !=0);
	
	
	#### caluate the average of deleterious mutation load and frequency of deleteriou mutants
	my $up_number=vector (@mut);
	my $up_average=median($up_number);
	my $up_mean=mean(@mut);
	
	### Also the standard deviation  
	my $stddev= stddev($up_number);
	### the max and min mutation load of deleterious mutation
	my $up_min=min(@mut);
	my $up_max=max(@mut);


	
	#### caluate the average of co-mutated mutation load and frequency of co-muatated mutants
	my $down_number=vector (@Othedelmut);
    my $down_average=median($down_number);
	my $down_mean=mean (@Othedelmut);
	my $down_max=max(@Othedelmut);
	my $down_min=min(@Othedelmut);
	
	
	
	

	
	#### caluate the average of co-mutated mutation load and frequency of co-muatated mutants
	my $Noncomut_number=vector (@Noncomut);
    my $Noncomut_average=median($Noncomut_number);
	my $Noncomut_mean=mean (@Noncomut);
	my $Noncomut_max=max(@Noncomut);
	my $Noncomut_min=min(@Noncomut);
	
	my $wilcoxon_test_p2= Statistics::Test::WilcoxonRankSum->new();
    my $ttest_p2=new Statistics::TTest;

	### put deleterious genes that only have one mutation across all the samples, it would be not possible to calcualte the pvalue
	
	my $pf1; my $pf1_ttest; my $pf2; my $pf2_ttest; 
	

	
	if ( $mu<2 || $Omut<2 || $Noncomut_num<2|| $up_average==0 || $down_average ==0 || $Noncomut_average==0){
		
		### gene 
		print OUT "$n\t-1\t-1\t-1\t-1\t$up_mean\t$down_mean\t$Noncomut_mean\t$mu\t$Omut\t$Noncomut_num\t$up_average\t$down_average\t$Noncomut_average\t$up_max\t$down_max\t$Noncomut_max\t$up_min\t$down_min\t$Noncomut_min\n";
		next;  ### this is to avoid the number is zero..
	}
	
	
	#### Step 4.1 calcualte P1
	#### Then measure the Pvalue of mutation loads from different samples
    my @mut_log=ln(@mut);
    my @non_mut_log=ln(@Othedelmut);

    $wilcoxon_test->load_data(\@mut_log, \@non_mut_log);
    $ttest->load_data(\@mut_log, \@non_mut_log);
	my $prob=$wilcoxon_test->probability();
	$pf1 = sprintf '%f', $prob;
	$pf1_ttest=sprintf '%f',$ttest-> {t_prob};
	
	
	

	
	##### Step 4.2 calcualte P2
	
    my @noncomut_log=ln(@Noncomut);
    $wilcoxon_test_p2->load_data(\@mut_log, \@noncomut_log);
    $ttest_p2->load_data(\@mut_log, \@noncomut_log);
	my $prob_2=$wilcoxon_test_p2->probability();
	$pf2 = sprintf '%f', $prob_2;
	$pf2_ttest=sprintf '%f',$ttest_p2-> {t_prob};
	
	
	
	
	print OUT "$n\t$pf1\t$pf2\t$pf1_ttest\t$pf2_ttest\t$up_mean\t$down_mean\t$Noncomut_mean\t$mu\t$Omut\t$Noncomut_num\t$up_average\t$down_average\t$Noncomut_average\t$up_max\t$down_max\t$Noncomut_max\t$up_min\t$down_min\t$Noncomut_min\n";
   	

}

close OUT;

sub mean {
	my $q=sprintf ('%f',sum(@_)/@_);
	return ($q);
}


sub ln{
        my @n=@_;
        my @t=();
        foreach my $n (@n){
                # if ($n==0){
                #       $n=0.01;
                # }
                my $t=sprintf ('%f', log($n+0.1)/log(10));
                push (@t,$t);
        }
        return (@t);
}

