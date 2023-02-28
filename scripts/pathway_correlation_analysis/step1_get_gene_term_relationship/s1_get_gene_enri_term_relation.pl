#!/usr/bin/perl -w
use warnings;

#usage: perl .pl map_id dir_name(IN_enricher_pathway_output) > gene_enri_term_relation.csv
#input file should only contain lf not crlf
($shared_id,$dir_name) = @ARGV ;
my$DIR_PATH = "./$dir_name"; #set the dir path here

opendir DIR, ${DIR_PATH} or die "can not open dir \"$DIR_PATH\"\n";
my@filelist = readdir DIR;
my@table = ();
my@table1 = ();
my$i = 0;

open IN, "$shared_id" or die;
while(<IN>){
	chomp;
	s/\r//g;
	s/"//g;
	$table[$i] = (split /\s+/,$_)[0];
	$table1[$i] = (split /\s+/,$_)[0];
	$i++;
}
#print "ID";

foreach my$file (@filelist) {
if($file =~ m/.csv/){
	my%hash = ();
	open II, "$DIR_PATH/$file" || die "can not open";
#	print "\t$`";
	while (<II>) {
		chomp;
		s/\r//g;
		s/"//g;
		my$id=(split /,/,$_)[0]; #<------target 
		my$data=(split /,/,$_)[1]; #<------target 
		$hash{$id}=$data;
	}
	close II;
	foreach(@table){
		chomp;
		my$id2=(split /,/,$_)[0];
		if(exists($hash{$id2})){
			$_ = "$_,\/$hash{$id2}";
		}else{
			$_ = "$_,\/0";
		}
	}
	}
}
#print "\n";
open OU, ">./tem1.txt";
foreach (@table)
{
	chomp;
#	s/\n/\t/;
	print OU "$_\n";
}
close OU;

open IN1, "./tem1.txt";
open OU1, ">./tem2.txt";
%hash1 = ();
while(<IN1>) {
	chomp;
	s/,/\t/;
	$id3 = (split /\t/,$_)[0];
	$data1 = (split /\t/,$_)[1];
	$data1 =~ s/,//g;
	$data1 =~ s/\/0//g;
	$data2 = $data1;
	$data1 =~ s/\/K/\nK/g;
	$hash1{$id3} = $data2;
	print OU1 "$data1";
}
print OU1 "\n";
close IN1;
close OU1;
#foreach $_(keys %hash1)
#{
#	print "$_\t$hash1{$_}\n";
#}

open IN2, "./tem2.txt";
%hash2 = ();
<IN2>;
while(<IN2>) {
	chomp;
	#print "$_\n";
	$hash2{$_} = 1;
}
close IN2;
#foreach $_(keys %hash2)
#{
#	print "$_\t$hash2{$_}\n";
#}

print "id";
foreach $_(keys %hash2)
{
	$id5 = $_;
	print ",$id5";
	foreach(@table1){
		chomp;
		my$id4=(split /,/,$_)[0];
		if($hash1{$id4} =~ m/$id5/){
			$_ = "$_,1";
		}else{
			$_ = "$_,0";
		}
	}
}
print "\n";
foreach (@table1)
{
	chomp;
	print "$_\n";
}
print "\n";

system("rm ./tem1.txt");
system("rm ./tem2.txt");




