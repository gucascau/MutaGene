#!/usr/bin/perl
#author:wangxin
### function: To determine the gene existence status between the low mutation samples and high mutation samples
### This scripts can also identify the driving genes that showed the greatly difference between two populations nonsilent and nonsilent

use strict;
use warnings;
use Statistics::Test::WilcoxonRankSum;
use List::Util qw(sum);

use Statistics::Multtest qw(bonferroni holm hommel hochberg BH BY qvalue);
use Statistics::Multtest qw(:all);
use Statistics::TTest;

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i}) {
       	die "************************************************************************
       	Usage: $0.pl
       		-i: query files (cancer type)
************************************************************************\n";
}



my $in=$opts{i};

open IN, "$in" or die $!;

while (<IN>){
	chomp;
	next if (/^#/);

   	my $cancer=$_;
	
	mkdir ($cancer);
	chdir ($cancer);
	
	
	system ("perl  /Users/xinwang/Documents/Project/Proposal/Figure3/scripts/Mutationload/Renew_P1P2_forcomuators.pl -i /Users/xinwang/Documents/Project/Proposal/Figure3/Version_20200511/All_mutationload_filtered_v0511.txt -g /Users/xinwang/Documents/Project/Proposal/Published_Data/mc3.v0.2.8.PUBLIC.maf -s $cancer");
	#system ("awk '{if(\$2<=0.01 && \$5>=10){print}}' $cancer.nonsil_all.txt > $cancer.nonsil_all_p0.01_n10.txt");
	
	
	chdir ("../");
}


close IN;



