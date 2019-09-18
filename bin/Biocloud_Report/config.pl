use Config::General;


my $config=&readConfig('detail.cfg');
&writeConfig($config,'aa.cfg');



sub writeConfig{
    my($config,$outFile)=@_;
#    if(-e $outFile){
#        system(qq(mv $outFile  $outFile.bak));
#    }
    my $con=Config::General->new(-SplitPolicy=> 'custom' , -SplitDelimiter => '\s*=\s*',-StoreDelimiter=>'  ',-SaveSorted => 1);
     $con->save_file($outFile, \%$config);
}
sub readConfig{
	my$configFile=shift;
    my $d=Config::General->new(-ConfigFile => "$configFile");
    my %config=$d->getall;
    return \%config;
}
