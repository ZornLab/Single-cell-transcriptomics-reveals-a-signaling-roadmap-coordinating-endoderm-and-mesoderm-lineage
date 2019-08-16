################## Script to create normalized MetaProfile
################## version 1
################## Script Requires two Inputs : 1) Directory with counts files [For example: Expression_Ligands_BMP.txt [this file has all BMP ligands and their counts] and similarly Expression_Receptors_Hedgehog.txt and so on] 2) OutputDirectory
################## How to run : perl GenerateMetaprofile_ForGeneSets.pl CountsDIR MetaProfileDIR
################## Please see that script requires counts files to be named in the following manner : Expression_Ligands_BMP.txt, Expression_Ligands_RA.txt, Expression_Response_FGF.txt etc.

use POSIX;
use List::Util qw(sum);
use List::Util qw( min max );
my $dir = $ARGV[0]; chomp $dir;                  #############----> Directory with Counts Files of GeneSets
my $resdir = $ARGV[1]; chomp $resdir;            #############----> OutPut Directory

if(-d $resdir)
{
 print "$resdir is already present \.\.\. Will use the same to write the results\n";
}
else
{
 print "Creating $resdir folder to write results\n";
 mkdir($resdir);
}

my $header_new;
open(OUT,">$resdir/MetaProfile_GeneSets.txt");
opendir(F,"$dir");
foreach my $f(readdir(F))
{
 if($f =~ /.txt/)
 {
  my @arrf=split(/\.txt/,$f);
  open(R,"$dir/$f");
  my $counter=0; my @header; my $numberOfSamples; my %hash_exp;
  while(my $data = <R>)
  {
   chomp $data;
   $counter++;
   if($counter==1)
   {
    $header_new=$data;
    break;
   }
  }
  close R;
 }
}
closedir F;
print OUT "$header_new\n";
opendir(F,"$dir");
foreach my $f(readdir(F))
{
 if($f =~ /.txt/)
 {
  my @arrf=split(/\.txt/,$f);
  my @arrf_newv=split(/\_/,$arrf[0]);
  my @arrf_newv1=split(//,$arrf_newv[1]);
  my $metaname=$arrf_newv[2]."\-".$arrf_newv1[0].$arrf_newv1[1].$arrf_newv1[2];
  print "Reading $f\n";
  open(R,"$dir/$f");
  my $counter=0; my @header; my $numberOfSamples; my %hash_exp;
  while(my $data = <R>)
  {
   chomp $data;
   $counter++;
   if($counter==1)
   {
    my @arr=split(/\t/,$data,2);
    @header=split(/\t/,$arr[1]);
    $numberOfSamples=scalar(@header);
   }
   else
   {
    my @arr=split(/\t/,$data,2);
    my @exp=split(/\t/,$arr[1]);
    my $max = abs(max(@exp));
    if($max == 0)
    {
     $hash_exp{$arr[0]}=$arr[1];
    }
    if($max > 0)
    {
     my @normalizedexp = map { $_/$max } @exp;
     my $normexp=join("\t",@normalizedexp);
     $hash_exp{$arr[0]}=$normexp;
    }
   }
  }
  close R;
  print "Strored all the information from $f\n";
  my $finalexp1="";
  for(my $i=0;$i<@header;$i++)
  {
   my $sum=0;
   foreach my $m(sort keys %hash_exp)
   {
    my @exp_get=split(/\t/,$hash_exp{$m});
    $sum=$sum+$exp_get[$i];
   }
   $finalexp1=$finalexp1."\t".$sum;
  }
  my @arr_new_exp=split(/\t/,$finalexp1);
  shift(@arr_new_exp);
  my $max_new=abs(max(@arr_new_exp));
  if($max_new == 0)
  {
   my $finalmetaexp=join("\t",@arr_new_exp);
   print OUT "$metaname\t$finalmetaexp\n";
   print "$metaname processed $f\n";
  }
  if($max_new > 0)
  {
   my @normalizedexpFinal = map { $_/$max_new } @arr_new_exp;
   my $finalmetaexp=join("\t",@normalizedexpFinal);
   print OUT "$metaname\t$finalmetaexp\n";
   print "$metaname processed $f\n";
  }
 }
}
closedir F;
close OUT;