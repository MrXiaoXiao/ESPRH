#!/usr/bin/perl -w
$year = "2017";
$mon = "02";
$day = "01";

$D = "$year/$mon/$day";
$R = "0.2/100/0.03/5/0.01";
$G = "0.5/100/0.01/5";
$V = "6.3/3.4/2.0/1.5/0";
$S = "3/1/4/1/2.0/0/5.0/5.0";

$dir = "./PickNetPicks/";
$station = "./sta_info_real_format.dat";
$ttime = "../REAL_scripts/tt_db/ttdb.txt";

system("REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime");
print"REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime\n";