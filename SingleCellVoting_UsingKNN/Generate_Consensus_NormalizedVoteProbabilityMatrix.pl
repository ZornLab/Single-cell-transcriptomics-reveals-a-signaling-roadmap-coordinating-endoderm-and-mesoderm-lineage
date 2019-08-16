############## Perl Script to generate Consensus Normalized vote probability matrix
############## version1
############## Script requires two inputs : 1) ProbabilityMatrix.txt [from KNN_Classification.r] 2) MetaFile.txt [Training Set Cells and their cluster annotation]
############## how to run : perl Script_Generate_Consensus_Matrix.pl ProbabilityMatrix.txt MetaFile.txt 


use List::Util 'max';
my $file = $ARGV[0]; chomp $file; 	########### ProbabilityMatrix.txt [from KNN_Classification.r]
my $file1 = $ARGV[1]; chomp $file1; ########### MetaFile.txt [Training Set Cells and their cluster annotation]
my $counternew=0; my %hash_traindata_clusters; my %hash_clustercells_train;
open(F,"$file1");
while(my $data = <F>)
{
 chomp $data;
 $counternew++;
 if($counternew==1)
 {}
 else
 {
  my @arr=split(/\t/,$data);
  $hash_traindata_clusters{$arr[1]}{$arr[0]}=0;
 }
}
close F;
foreach my $m(sort keys %hash_traindata_clusters)
{
 my $countcells_train=0;
 foreach my $n(keys %{$hash_traindata_clusters{$m}})
 {
  $countcells_train++;
 }
 $hash_clustercells_train{$m}=$countcells_train;
}
my $counter=0; my @header; my %hash_store; my %hash_res; my $header_n;
open(F,"$file");
while(my $data = <F>)
{
 chomp $data;
 $counter++;
 if($counter==1)
 {
  my @header_temp=split(/\t/,$data,3);
  my @header_temp1=split(/\t/,$header_temp[2]);
  foreach my $m(@header_temp1)
  {
   if($m eq "")
   {}
   else
   {
    push(@header,$m);
   }
  }
 }
 else
 {
  my @arr=split(/\t/,$data,3);
  my @arr1=split(/\t/,$arr[2]);
  for(my $i=0;$i<@arr1;$i++)
  {
    my $value=$arr[0]."\t".$arr1[$i];
    my $key=$arr[1]."\t".$header[$i];
    $hash_store{$key}{$value}=0;
  }
 }
}
close F;
foreach my $m(sort keys %hash_store)
{
 my @arr=split(/\t/,$m);
 my $sum_of_votes=0; my $count_cells=0;
 foreach my $n(keys %{$hash_store{$m}})
 {
  my @arr1=split(/\t/,$n);
  $count_cells++;
  $sum_of_votes=$sum_of_votes+$arr1[1];
 }
 my $normalizedProb=$sum_of_votes/$count_cells;
 if(exists $hash_clustercells_train{$arr[1]})
 {
  my $valuetemp=($normalizedProb/$hash_clustercells_train{$arr[1]})*1000; ########## Changes done
  my $value=$arr[1]."\t".$valuetemp;
  $hash_res{$arr[0]}{$value}=0;
 }
}
foreach my $m(@header)
{
 $header_n=$header_n."\t".$m;
}
open(OUT,">NormalizedConsensus_VoteProbability.txt");
print OUT "Cluster$header_n\n";
foreach my $m(sort keys %hash_res)
{
 my $final_value="";
 foreach my $a(@header)
 {
  my $val="";
  foreach my $n(sort keys %{$hash_res{$m}})
  {
   my @arr=split(/\t/,$n);
   if($a eq $arr[0])
   {
    $val=$arr[1];
   }
  }
  if($val eq "")
  {
   $final_value=$final_value."\t"."NA";
  }
  else
  {
   $final_value=$final_value."\t".$val;
  }
 }
 print OUT "$m$final_value\n";
}
close OUT;