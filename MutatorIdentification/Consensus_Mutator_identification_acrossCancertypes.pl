#!/usr/bin/perl
#author:wangxin
### function: Here we also restricted the sample of specific gene mutations to silent only. For the same candidate genes,  we calculate the significance between silent-only mutatants and no-mutants.

use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i}|| !defined $opts{o} ) {
       	die "************************************************************************
       	Usage: $0.pl
	       		-i: query of regulator folder including cancer types
			-o: output with the information of existence status in each cancer, another with average mutation load in Nonsil, Sil and Nomutation
************************************************************************\n";
}


my $input=$opts{i};
my $out=$opts{o};


my @dir;
opendir IN,"$input" or die $!;
	@dir=readdir(IN);
closedir IN;

my %hash; my %Sil; my %NonSil; my %NonMut; my %mut;

foreach my $i (@dir){
	if ($i=~/(\S+)\.mutationload.p1p2.txt.mutators.txt/){
		my $type=$1;
		open I,"$i" or die $!;
		while (<I>){
			chomp;
			my @array=split/\t/,$_;
			$hash{$array[0]}->{$type}++;
			$mut{$type}++;
# 			$Sil{$array[0]}->{$type}=$array[2];
# 			$NonSil{$array[0]}->{$type}=$array[8];
# 			$NonMut{$array[0]}->{$type}=$array[3];
		}
		close I;
	}
	
}


open OUT, ">$out.exist" or die $!;
#open M,">$out.mutation" or die $!;
print OUT "Gene";
#print M "Gene";
foreach my $i (sort keys %mut){
	print OUT "\t$i";
	#print M "\t$i.Nonsil\t$i.Sil\t$i.NonMut";
}
print OUT "\n";
#print M "\n";


foreach my $g (sort keys %hash){
	print OUT "$g";
	#print M "$g";
	foreach my $t (sort keys %mut){
		$hash{$g}->{$t}=($hash{$g}->{$t})?1:0;
		# $Sil{$g}->{$t}=($Sil{$g}->{$t})?$Sil{$g}->{$t}:0;
# 		$NonSil{$g}->{$t} = ($NonSil{$g}->{$t})?$NonSil{$g}->{$t}:0;
# 		$NonMut{$g}->{$t} = ($NonMut{$g}->{$t})?$NonMut{$g}->{$t}:0;
#
		print OUT "\t$hash{$g}->{$t}";
		#print M "\t$NonSil{$g}->{$t}\t$Sil{$g}->{$t}\t$NonMut{$g}->{$t}";
	}
	print OUT "\n";
	#print M "\n";
}

close OUT;

#close M;
