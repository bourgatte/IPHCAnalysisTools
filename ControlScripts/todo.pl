#! /usr/bin/perl
use Cwd;
use POSIX;
use POSIX qw(strftime);

#############################################
$numArgs = $#ARGV +1;
$ARGV[$argnum];

$UserID= POSIX::cuserid();
$UserIDCern=$UserID;
$UserDir="";
$UserName="";
if($UserID eq "vcherepa"){
    $UserIDCern="cherepan";
    $UserDir="--cherepanov";
    $UserName="Vladimir";
}

if($UserID eq "goe"){
    $UserIDCern="goe";
    $UserDir="--goe";
    $UserName="Ulrich";
}
if($UserID eq "lebihan"){
    $UserIDCern="lebihan";
    $UserDir="--lebihan";
    $UserName="Anne-Catherine";
}

if($UserID eq "gbourgat"){
    $UserIDCern="gbourgat";
    $UserDir="--gbourgat";
    $UserName="Guillaume";
}


if($UserID eq "cgrimault"){
    $UserIDCern="cgrimault";
    $UserDir="--cgrimault";
    $UserName="Clement";
}

if($UserID eq "msessini"){
    $UserIDCern="msessini";
    $UserDir="--msessini";
    $UserName="Mario";
}

#Default values
$InputDir="/home-pbs/$UserID/InputTest";
$OutputDir="/home-pbs/$UserID/Test";
#$OutputDir="~/";
$CodeDir="../Code";
$set="ControlSample_";
$CMSSWRel="8_0_26_patch1";
$PileupVersion="V08-03-17";
$tag="03-00-12";
$TauReco="5_2_3_patch3_Dec_08_2012";
$BTag="NO";
$Cleaning ="NO";
$maxdata=20;
$maxmc=5;
$maxemb=20;
$ARCH="slc6_amd64_gcc530";
$Queue="cream-pbs-short";
$QsubQue="cms_local_short";


if($ARGV[0] eq "--help" || $ARGV[0] eq ""){
    printf("\n ===========>>>>  Hello  $UserName ! <<<<=============\n\n");
    printf("\nThis code requires one input option. The systax is:./todo_Grid.pl [OPTION]");
    printf("\nPlease choose from the following options:\n");
    printf("\n./todo.pl --help                                   Prints this message\n");
    printf("\n./todo.pl --TauNtuple <TauNtupleDir>   --CMSSWRel <X_Y_Z>            This Option will perform installation of LLR Ntuple framework production");
    printf("\n./todo.pl --Local <Input.txt>                      INTENTED FOR SMALL SCALE TESTS ONLY");  
    printf("\n                                                   Configure a directory to run locally. <InputPar.txt> name of file that");
    printf("\n                                                   contains input command template.");
    printf("\n                                                   Optional commands:  ");
    printf("\n                                                     --InputDir <InputDir>   Default value: $InputDir"); 
    printf("\n                                                       Note: the root files must be inside a $InputDir/<dir>/ ");
    printf("\n                                                       where dir contains user, data or mc as a substring");
    printf("\n                                                     --OutputDir <OutputDir> Default value: $OutputDir");
    printf("\n                                                     --CodeDir  <CodeDir>    Default value: $CodeDir");
    printf("\n                                                     --SetName <SetName>     Default value: $set ");
    printf("\n                                                     --NMaxData <Max Number of data files per job >     Default value: $maxdata ");
    printf("\n                                                     --NMaxMC <Max Number of MC files per job >     Default value: $maxmc ");
    printf("\n                                                     --NMaxEmbed <Max Number of Embedding files per job > Default value: $maxemb ");
    printf("\n                                                     --ROOTSYS <ROOTSYS> the current ROOTSYS variable if --BuildRoot is not defined");
    printf("\n                                                     --TauSpinner Option to turn on TauSpinner");
    printf("\n                                                     --SVfit Option to turn on SVfit");
    printf("\n                                                     --QsubQueue <queue>; type of queue of jobs submitted by qsub, example:  ");
    printf("\n                                                     --QsubQueue short, --QsubQueue medium, --QsubQueue long; Default: short.");
    printf("\n                                                     --Proxy <path>;   a path to proxy recieved by running grid-init");
    printf("\n./todo.pl --DCache <Input.txt> <ListofDS.txt>      INTENTED FOR REGULAR USE (DEFAULT)");
    printf("\n                                                   Configure a directory to run from. <InputPar.txt> name of file that");
    printf("\n                                                   contains input command template.");
    printf("\n                                                   <ListoDS.txt> list of DCache Dataset directories you want to run on.");
    printf("\n                                                   Optional commands:  ");
    printf("\n                                                     --OutputDir <OutputDir> ");
    printf("\n                                                     --CodeDir <CodeDir>");
    printf("\n                                                     --SetName <SetName> "); 
    printf("\n                                                     --NMaxData <Max Number of data files per job >     Default value: $maxdata ");
    printf("\n                                                     --NMaxMC <Max Number of MC files per job >     Default value: $maxmc ");
    printf("\n                                                     --NMaxEmbed <Max Number of Embedding files per job > Default value: $maxemb ");
    printf("\n                                                     --ROOTSYS <ROOTSYS> the current ROOTSYS variable if --BuildRoot is not defined");
    printf("\n                                                     --TauSpinner Option to turn on TauSpinner");
    printf("\n                                                     --SVfit Option to turn on SVfit --- Not implemented at the moment");
    printf("\n                                                     --QsubQueue <queue>; type of queue of jobs submitted by qsub, example:  ");
    printf("\n                                                     --QsubQueue short, --QsubQueue medium, --QsubQueue long; Default: short.");
    printf("\n                                                     --Proxy <path>;   a path to proxy recieved by running grid-init");
    printf("\n  ");
    printf("\n./todo.pl --GRID <Input.txt> <ListofDS.txt> --ROOTSYS \$ROOTSYS       ALTERNATIVE FOR REGULAR USE  (will be implemented later)");
    printf("\n                                                   Configure a directory to run from. <InputPar.txt> name of file that");
    printf("\n                                                   contains input command template.");
    printf("\n                                                   <ListoDS.txt> list of DCache Dataset directories you want to run on.");
    printf("\n                                                     --ROOTSYS <ROOTSYS> the current ROOTSYS variable if --BuildRoot is not defined");
    printf("\n                                                   Optional commands:  ");
    printf("\n                                                     --OutputDir <OutputDir> ");
    printf("\n                                                     --CodeDir <CodeDir>");
    printf("\n                                                     --SetName <SetName> ");
    printf("\n                                                     --NMaxData <Max Number of data files per job > Default value: $maxdata");
    printf("\n                                                     --NMaxMC <Max Number of MC files per job > Default value: $maxmc ");
    printf("\n                                                     --NMaxEmbed <Max Number of Embedding files per job > Default value: $maxemb ");
    printf("\n                                                     --BuildRoot <ROOT Version> builds custom version for root instead of copying lib+include");
    printf("\n                                                     --GRIDSite <site> the grid site you wish to run on. Default=grid-srm.physik.rwth-aachen.de");
    printf("\n                                                     --LongQueue Option to run on CMS queue (allows for longer jobs)");
    printf("\n                                                     --TauSpinner Option to turn on TauSpinner");
    printf("\n                                                     --SVfit Option to turn on SVfit \n ");
    exit(0);  
} 

######################################
$InputFile=$ARGV[1];
$buildRoot=0;
$hasroot=0;
#$ubergridsite="grid-ftp.physik.rwth-aachen.de";
#$dcapgridsite="grid-dcap.physik.rwth-aachen.de";
#$gridsite="grid-srm.physik.rwth-aachen.de";
$tauspinner="";
$svfit="";
for($l=2;$l<$numArgs; $l++){
    if($ARGV[$l] eq "--InputDir"){
	$l++;
	$InputDir=$ARGV[$l];
    }
    if($ARGV[$l] eq "--OutputDir"){
	$l++;
	$OutputDir=$ARGV[$l];
    }
    if($ARGV[$l] eq "--CodeDir"){
	$l++;
	$CodeDir=$ARGV[$l];
    }
    if($ARGV[$l] eq "--SetName"){
	$l++;
	$set=$ARGV[$l];
    }
    if($ARGV[$l] eq "--CMSSWRel"){
        $l++;
        $CMSSWRel=$ARGV[$l];
    }
    if($ARGV[$l] eq "--PileupVersion"){
       $l++;
        $PileupVersion=$ARGV[$l];
    }
    if($ARGV[$l] eq "--tag"){
        $l++;
        $tag=$ARGV[$l];
    }
    if($ARGV[$l] eq "--TauReco"){
        $l++;
        $TauReco=$ARGV[$l];
    }
    if($ARGV[$l] eq "--BTag"){
	$l++;
	$BTag=$ARGV[$l];
    }
    if($ARGV[$l] eq "--Cleaning"){
	$l++;
	$Cleaning =$ARGV[$l];
    }
    if($ARGV[$l] eq  "--NMaxData" ){
	$l++;
	$maxdata=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--NMaxMC" ){
	$l++;
	$maxmc=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--NMaxEmbed" ){
	$l++;
	$maxemb=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--BuildRoot"){
	$l++;
	$buildrootversion=$ARGV[$l];
	$buildRoot=1;
    } 
    if($ARGV[$l] eq  "--ROOTSYS"){
        $l++;
        $MYROOTSYS=$ARGV[$l];
	$hasroot=1;
    }
    if($ARGV[$l] eq  "--GRIDSite"){
	$l++;
	$ubergridsite="grid-ftp."+$ARGV[$l];
	$dcapgridsite="grid-dcap."+$ARGV[$l];
	$gridsite="grid-srm."+$ARGV[$l];
    }
    if($ARGV[$l] eq  "--LongQueue"){
		$l++;
		$Queue="cream-pbs-cms";
    }
    if($ARGV[$l] eq  "--ARCH" ){
        $l++;
        $ARCH=$ARGV[$l];
    }
   if($ARGV[$l] eq  "--TauSpinner" ){
	$tauspinner="--TauSpinner";
    }
    if($ARGV[$l] eq  "--Proxy" ){
        $l++;
        $Proxy=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--SVfit" ){
    $svfit="--SVfit";
    }
    if($ARGV[$l] eq  "--QsubQueue" ){
	$l++;
	if($ARGV[$l] eq  "short"){
	    $QsubQue="cms_local_short";
	}
	if($ARGV[$l] eq  "medium"){
	    $QsubQue="sbg_local_mdm";
	}
	if($ARGV[$l] eq  "long"){
	    $QsubQue="cms_local";
	}
    }
}
my $dir = getcwd;

$time= strftime("%h_%d_%Y",localtime);
$temp= $set . $time;
$set=$temp;

if( $ARGV[0] eq "--TauNtuple"){
    $currentdir=getcwd;
    if($ARGV[1] ne ""){
	$basedir=$ARGV[1];
    }
    else{
	printf("\nWorkingDir for CMSSW is required. Please follow the syntax:./todo.pl --TauNtuple <TauNtupleDir> ");
	printf("\nFor more details use: ./todo --help\n"); 
	exit(0);
    }
    printf("\nWorkingDir for CMSSW: $basedir");
    printf("\nCurrentDir is: $currentdir \n");

    system(sprintf("rm Install_TauNtuple_$time"));

    system(sprintf("echo \"cernlib-use root\" >> Install_TauNtuple_$time"));

    system(sprintf("echo \"export SCRAM_ARCH=\\\"$ARCH\\\"\" >> Install_TauNtuple_$time"));
    system(sprintf("echo \"source /cvmfs/cms.cern.ch/cmsset_default.sh\" >> Install_TauNtuple_$time"));


    system(sprintf("echo \"cd $basedir\" >>  Install_TauNtuple_$time")); 
    system(sprintf("echo \"cmsrel CMSSW_$CMSSWRel\" >>  Install_TauNtuple_$time")); 
    system(sprintf("echo \"cd CMSSW_$CMSSWRel/src\" >> Install_TauNtuple_$time")); 
    system(sprintf("echo \"cmsenv\" >> Install_TauNtuple_$time")); 
    $CMSPATH="/CMSSW_$CMSSWRel/";
    $CMSSW_BASE="$basedir$CMSPATH";
   system(sprintf("echo \"git cms-merge-topic cms-met:METRecipe_8020\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git cms-merge-topic ikrav:egm_id_80X_v2\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git cms-merge-topic gpetruc:badMuonFilters_80X_v2\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone https://github.com/cherepan/LLRHiggsTauTau\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd LLRHiggsTauTau; git checkout VladimirDev\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd -\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd -\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd -\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd UFHZZAnalysisRun2 ; git checkout master FSRPhotons\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd -\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd TauAnalysis/SVfitStandalone\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git checkout svFit_2015Apr03\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd $CMSSW_BASE/src\" >> Install_TauNtuple_$time")); 

   system(sprintf("echo \"scram b -j 16\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd $CMSSW_BASE/external/$SCRAM_ARCH\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd data/RecoEgamma/ElectronIdentification/data\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git checkout egm_id_80X_v1\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd $CMSSW_BASE/src\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"scram b -j 16\" >> Install_TauNtuple_$time")); 
    # print Instructions
    printf("\n\nInstructions:");
    printf("\nsource  Install_TauNtuple_$time to complete installation, compilation might take some time...  \n\n");
    printf("\nTo run test job do  'cmsRun analyzer.py'  in  $CMSSW_BASE/src/LLRHiggsTauTau/NtupleProducer/test \n\n");

}
 
if( $ARGV[0] eq "--Local" ){
    # Print out input parameters
    printf("Active directory will be: $OutputDir/workdir$set \n");
    printf("Input directory will be:  $InputDir \n");
    printf("Code Repository is:       $CodeDir \n");
    # Clean Directory in case workdir$set exists
    printf("Cleaning Directories \n");
    system(sprintf("cd $OutputDir"));
    system(sprintf("rm -rf $OutputDir/workdir$set \n"));
    system(sprintf("mkdir $OutputDir/workdir$set "));
    printf("Cleaning complete \n");
    
    # create directory stucture
    system(sprintf("mkdir  $OutputDir/workdir$set/Code "));
    system(sprintf("mkdir  $OutputDir/workdir$set/Code/i386_linux "));
    system(sprintf("cp -r $CodeDir/* $OutputDir/workdir$set/Code/ "));
    system(sprintf("mkdir $OutputDir/workdir$set/EPS "));
    system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/InputData "));
    
    # generate compile script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"source config \\\$\@ \"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake all  \" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/compile")) ;
 
    # Generate Combine script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Combine")) ;
    system(sprintf("echo \"export workdir=\\\"$OutputDir/workdir$set/\\\"\" >> $OutputDir/workdir$set/Combine"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Combine")); 
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Combine")) ; 
    system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe \" >> $OutputDir/workdir$set/Combine")) ;

    # Generate Combine Input
    system(sprintf("cp $InputFile $OutputDir/workdir$set/Input.txt "));
    system(sprintf("cd $OutputDir/workdir$set; $dir/subs '{SET}' COMBINE Input.txt; cd $dir "));
    system(sprintf("cd $OutputDir/workdir$set; $dir/subs '{FileDir}' NotAvailable Input.txt; cd $dir "));
    system(sprintf("echo \"Mode: RECONSTRUCT\" >> $OutputDir/workdir$set/Input.txt"));
    system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Input.txt"));

    # Setup Condor Combine scripts
#    system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Condor_Combine"));
#    system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Condor_Combine"));
#    system(sprintf("echo \"executable   = Combine      \"  >> $OutputDir/workdir$set/Condor_Combine")); 
#    system(sprintf("echo \"output       = Combine-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Condor_Combine")); 
#    system(sprintf("echo \"error        = Combine-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Condor_Combine")); 
#    system(sprintf("echo \"log          = Combine-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Condor_Combine")); 			
#    system(sprintf("echo \"queue = 1 \" >> $OutputDir/workdir$set/Condor_Combine"));

    # Start Submit script
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Submit")) ; 
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"rm Set*/*.o; rm Set*/*.e; rm Set*/*.log; \" >> $OutputDir/workdir$set/Submit")) ;

    opendir(DIR,"$InputDir/");
    @dirs = grep {( /user/ || /data/ || /mc/)} readdir(DIR);
    closedir DIR;
    
    $B=0;
    $max=1;
    for($l=0;$l<2; $l++){
	printf("\n\nStarting Loop $l  maxdata = $maxdata maxmc = $maxmc maxemb = $maxemb \n");
	$A=$maxdata+$maxmc+$maxemb+10;
	foreach $subdir (@dirs){
	    printf("subdir  $subdir\n");
	    if(($l==0 && ($subdir =~ m/data/)) || ($l==1 && !($subdir =~ m/data/))){
		if($l==0){
		    $max=$maxdata;
		}
		else{
		    $max=$maxmc;
			if($subdir =~ m/embed/){
			    $max=$maxemb;
			}
		}
		printf(" \nAccessing Directory  $subdir \n");
		opendir(SUBDIR,"$InputDir/$subdir/");
		@files = grep { /root/ } readdir(SUBDIR);

		$nfiles = @files;
		$idx=0;
		foreach $file (@files){
		    $idx++;
		    printf("$file Set= $B  Index=  $A   Max.= $max N Files = $nfiles Current File = $idx \n");
		    if($A > $max){
			$A=1;
			$B++;
			# Add Set information to Combining scripts and Input.txt
			system(sprintf("echo \"File: $OutputDir/workdir$set/Set_$B/  \" >>  $OutputDir/workdir$set/Input.txt ")) ;
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B \" >> $OutputDir/workdir$set/Submit")) ;
			system(sprintf("echo \"source Qsub_Set_$B  \" >> $OutputDir/workdir$set/Submit")) ;
#			system(sprintf("echo \"condor_submit  Condor_Set_$B  \" >> $OutputDir/workdir$set/Submit")) ;

			# Create and configure Set_$B dir
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B ")) ;
			system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/Set_$B/InputData "));
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B/EPS ")) ;

                        # Setup Set_$B.sh
			system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("echo \"echo 'Starting Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"export workdir=\\\"$OutputDir/workdir$set/\\\"\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("chmod +x $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe 2>&1 | tee >(sed -r \\\"s/\\\\x1B\\\\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g\\\" > Set_$B.output) \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"echo 'Completed Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));

			# Setup Input.txt
			system(sprintf("cp   $InputFile $OutputDir/workdir$set/Set_$B/Input.txt ")); 
			system(sprintf("cd $OutputDir/workdir$set/Set_$B; $dir/subs '{SET}' Set_$B Input.txt; cd $dir "));
			system(sprintf("cd $OutputDir/workdir$set/Set_$B; $dir/subs '{FileDir}' $InputDir/$subdir Input.txt; cd $dir ")); 
			system(sprintf("echo \"Mode: ANALYSIS\" >> $OutputDir/workdir$set/Set_$B/Input.txt")); 
			system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Set_$B/Input.txt"));



			# Setup qsub scripts
			$s1_char='\${que}';
			$s2_char='\${output}';
			$s3_char='\${error}';
			system(sprintf("echo \" #PBS -u $UserID\"  >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B"));
			system(sprintf("echo \" #! /bin/bash\"  >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B"));
			system(sprintf("echo \" export HOME=\\\"$OutputDir/workdir$set/\\\"         \"  >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B")); 
			system(sprintf("echo \" que=\\\"$QsubQue\\\"\" >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B")); 
			system(sprintf("echo \" output=\\\"Set_$B.qsub.o  \\\"\" >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B")); 
			system(sprintf("echo \" error=\\\"Set_$B.qsub.e  \\\"\" >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B")); 
			system(sprintf("echo \" qsub -q  $s1_char  -o $s2_char -e $s3_char Set_$B.sh\" >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B"));		
		    }
		    system(sprintf("echo \"File: $InputDir/$subdir/$file  \" >> $OutputDir/workdir$set/Set_$B/Input.txt")) ;
		    $A++;
		}
	    }
	}
    }
    
    printf("Setting up root....\n");
    if($buildRoot==1){
        printf("Building custom root $buildrootversion... Enjoy your coffee!!! ");
	system(sprintf("wget ftp://root.cern.ch/root/root_v$buildrootversion.source.tar.gz"));
	system(sprintf("gzip -dc root_v$buildrootversion.source.tar.gz | tar -xf -"));
	system(sprintf("mkdir $OutputDir/workdir$set/root"));
	system(sprintf("cd root_v$buildrootversion; ./configure --enable-python --enable-roofit --enable-minuit2 --disable-xrootd --disable-sqlite --disable-python --disable-mysql --prefix=$OutputDir/workdir$set/root; make & make install "));
    }

    # Finish Submit script
    system(sprintf("echo \"cd  $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit"));

    # print Instructions  
    printf("\n\nInstructions");
    printf("\nPlease make sure you have run:");
    printf("\ngit config --global credential.helper 'cache --timeout=3600'");
    printf("\ngit config --global credential.helper cache");
    printf("\nNow you can run the analysis using dcache.");
    printf("\nTo go to the Test workdir: cd  $OutputDir/workdir$set ");
#    printf("\nTo compile the code in the workdir: source compile  --useRoot $OutputDir/workdir$set/root/ $UserDir $tauspinner $svfit");
    printf("\nTo compile the code in the workdir: source compile   $UserDir $tauspinner $svfit");
    printf("\nTo submit jobs to the batch queue: source Submit ");
    printf("\nTo combine jobs submitted to the batch queue: source Combine \n");
#    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile  --useRoot $OutputDir/workdir$set/root/ $UserDir; cd $OutputDir/workdir$set/Set_1; source Set_1 | tee log; cd ..\n");
    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile   $UserDir $tauspinner $svfit; cd $OutputDir/workdir$set/Set_1; source Set_1 | tee log; cd ..\n");

}

#xrdcp root://sbgse1.in2p3.fr//store/user/cherepan/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DefaultPublishName/170523_140637/0000/HTauTauAnalysis_425.root .


 if( $ARGV[0] eq "--DCache" ){
    $RemoteScrathDir="/opt/sbg/scratch1/cms/";
    $RemoteDir='\$TMPDIR'; 
    #$RemoteDir   <-> $RemoteDir$UserID
    $TempDataSetFile=$ARGV[2];
    # Print out input parameters
    printf("Active directory will be: $OutputDir/workdir$set \n");
    printf("Code Repository is:       $CodeDir \n");
    printf("List of dcache dir:       $TempDataSetFile \n");

    # Open ListofFile.txt
    @DataSets;
    open(DAT, $TempDataSetFile) || die("Could not open file $TempDataSetFile!");
    while ($item = <DAT>) {
	chomp($item);
	push(@DataSets,$item);
    }
    close(DAT);
    # Clean Directory in case workdir$set exists
    printf("Cleaning Directories \n");
    system(sprintf("cd $OutputDir"));
    system(sprintf("rm -rf $OutputDir/workdir$set \n"));
    system(sprintf("mkdir $OutputDir/workdir$set "));
    printf("Cleaning complete \n");
    
    # creat directory stucture
    system(sprintf("mkdir  $OutputDir/workdir$set/Code "));
    system(sprintf("mkdir  $OutputDir/workdir$set/Code/i386_linux "));
    system(sprintf("cp -r $CodeDir/* $OutputDir/workdir$set/Code/ "));
    system(sprintf("mkdir $OutputDir/workdir$set/EPS "));
    system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/InputData "));
    
    # generate compile script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"source config \\\$\@ \"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake all\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/compile")) ;
 
    # Generate Combine script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Combine")) ;
    system(sprintf("echo \"export workdir=\\\"$OutputDir/workdir$set/\\\"\" >> $OutputDir/workdir$set/Combine"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Combine"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Combine")) ; 
    system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe \" >> $OutputDir/workdir$set/Combine")) ;

    # Generate Combine Input
    system(sprintf("cp $InputFile $OutputDir/workdir$set/Input.txt "));
    system(sprintf("cd $OutputDir/workdir$set/; $dir/subs '{SET}' COMBINE Input.txt; cd $dir"));
    system(sprintf("cd $OutputDir/workdir$set/; $dir/subs '{FileDir}' COMBINE Input.txt; cd $dir"));
    system(sprintf("echo \"Mode: RECONSTRUCT\" >> $OutputDir/workdir$set/Input.txt"));
    system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Input.txt"));

    # Generate SubmitManual and set_env 
    system(sprintf("cp  subs  $OutputDir/workdir$set/;"));
    system(sprintf("cp  SubmitManual  $OutputDir/workdir$set/;"));
    system(sprintf("cd $OutputDir/workdir$set/; ./subs '{DIR}'  $OutputDir/workdir$set/  $OutputDir/workdir$set/SubmitManual; "));
			
    system(sprintf("cp  set_env  $OutputDir/workdir$set/;"));
    system(sprintf("cd $OutputDir/workdir$set/; ./subs '{DIR}'  $RemoteScrathDir$UserID/  $OutputDir/workdir$set/set_env; "));



    # Setup Condor Combine scripts
    #system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Condor_Combine"));
    #system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Condor_Combine"));
    #system(sprintf("echo \"executable   = Combine      \"  >> $OutputDir/workdir$set/Condor_Combine")); 
    #system(sprintf("echo \"output       = Combine-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    #system(sprintf("echo \"error        = Combine-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    #system(sprintf("echo \"log          = Combine-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Condor_Combine")); 			
    #system(sprintf("echo \"queue = 1 \" >> $OutputDir/workdir$set/Condor_Combine"));

    # Start Submit script
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Submit")) ; 
    system(sprintf("echo \"verbosity=\\\$(grep SetLevel Code/Analysis.cxx | grep -c -e Debug -e Verbose)\" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"  if [[ \\\${verbosity} -ne 0 ]]; then \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"    echo 'ERROR: Please make sure to set the verbosity level to Info in Analysis.cxx, otherwise your log-files will break QSUB! Abort...' \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"    exit \\\${verbosity}\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  fi\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"rm Set*/*.o; rm Set*/*.e; rm Set*/*.log; \" >> $OutputDir/workdir$set/Submit")) ;
 

    $B=0;
    for($l=0;$l<2; $l++){
	system(sprintf("echo \"notification = Error        \" >> $OutputDir/workdir$set/Condor_Combine"));
	$max=1;
	foreach $DS (@DataSets){
	    #if((($l==0 && ($DS =~ m/Data/)) || ($l==1 && !($DS =~ m/Data/))) || (($l==0 && ($DS =~ m/data/)) || ($l==1 && !($DS =~ m/data/)))){
	    if((($l==0 && ($DS =~ m/Data/)) || ($l==1 && !($DS =~ m/Data/)))){
		#	print "true 1   l = $l\n";
		if($l==0){
		    #	print "true 2   l = $l\n";
		    $max=$maxdata;
		}
		else{
		    #	print "true 3   l = $l\n";
		    $max=$maxmc;
		    if($DS =~ m/embed/){
			#	print "true 4";
			$max=$maxemb;
		    }
		}
		
		print "max  = $max;    maxmc = $maxmc \n";
#		print "------- DS   $DS  \n";
		printf("\n\nStarting Loop $l \n");
		$A=$maxdata+$maxmc+$maxemb+10;

		# find the root files for the current DataSet (DS)
		printf("Accessing Directory  $DS \n");
		system(sprintf("touch junk0"));
		system(sprintf("touch junk1"));
		system(sprintf("touch junk2"));
		system(sprintf("touch junk"));
		system(sprintf("rfdir /dpm/in2p3.fr/home/cms/phedex$DS >& junk0 "));
		@dpmsubdirs=();
		open(DAT, "junk0");
		while ($item = <DAT>) {
		    chomp($item);
		    push(@dpmsubdirs,$item);
		}
		close(DAT);
		system(sprintf("cat junk0 | awk '{print \$9}' >& junk1")); 
		@dpmdirs=();
		open(DAT, "junk1");
		while ($item = <DAT>) {
		    chomp($item);
		    $fpath="$DS$item";
		    push(@dpmdirs,$fpath);
		}
	#	print @dpmdirs;
		# Get list of files in dcache dir
		@files=(); 
		foreach $ipath (@dpmdirs){
		    printf(" path =  $ipath  \n ");
		    system(sprintf("rfdir /dpm/in2p3.fr/home/cms/phedex$ipath | grep \".root\" >& junk2 "));
		    system(sprintf("cat junk2 | awk '{print \$9}' >& junk")); 
		    open(DAT, "junk");
		    while ($item = <DAT>) {
			chomp($item);
#			printf(" file --- =  $ipath/$item  \n ");
			$FileList="$ipath/$item";
			#print ">>>>>>>>>>>> Filelst  $FileList";
			push(@files,$FileList);
		    }
		}

		close(DAT);
		system(sprintf("rm junk"));
		system(sprintf("rm junk0"));
		system(sprintf("rm junk1"));
		system(sprintf("rm junk2"));

		$nfiles = @files;
		$idx=0;
	    
		foreach $file (@files){
		    $idx++;
		    printf("$file Set= $B  Index=  $A   Max.= $max N Files = $nfiles Current File = $idx \n");
		    if($A > $max ){
			$A=1;
			$B++;

			# Add Set information to Combining scripts and Input.txt
			system(sprintf("echo \"File: $OutputDir/workdir$set/Set_$B/ \" >>  $OutputDir/workdir$set/Input.txt ")) ;
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B \" >> $OutputDir/workdir$set/Submit")) ;
			system(sprintf("echo \"source  Qsub_Set_$B  \" >> $OutputDir/workdir$set/Submit")) ;


			# Create and configure Set_$B dir
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B ")) ;
			system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/Set_$B/InputData "));
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B/EPS ")) ;
			
			# Setup Set_$B.sh
			system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("echo \"echo 'Starting Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"export workdir=\\\"$OutputDir/workdir$set/\\\"\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"export X509_USER_PROXY=\\\"$Proxy\\\"\"  >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"source $OutputDir/workdir$set/Set_$B/Set_$B-get.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ; 
			system(sprintf("chmod +x $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"mkdir $RemoteDir/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cp -r *    $RemoteDir/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd  $RemoteDir/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe 2>&1 | tee >(sed -r \\\"s/\\\\x1B\\\\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g\\\" > Set_$B.output) \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cp -r *  $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"source $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"rm -r   $RemoteDir/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));	
			system(sprintf("echo \"export HOME=\\\"/home/$UserID\\\"         \"   >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"echo 'Completed Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));


                        # Setup Set_$B_get.sh and Set_$B_clean.sh
			system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			system(sprintf("echo \"mkdir $RemoteDir \" >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			system(sprintf("echo \"cd $RemoteDir \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));

			system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));
			system(sprintf("echo \"cd $RemoteDir \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh")); 

			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			# Setup Input.txt
			system(sprintf("cp   $InputFile $OutputDir/workdir$set/Set_$B/Input.txt "));
			system(sprintf("cd $OutputDir/workdir$set/Set_$B; $dir/subs '{SET}' Set_$B Input.txt; cd $dir "));
			system(sprintf("cd $OutputDir/workdir$set/Set_$B; $dir/subs '{FileDir}' $DS Input.txt; cd $dir "));
			system(sprintf("echo \"Mode: ANALYSIS\" >> $OutputDir/workdir$set/Set_$B/Input.txt")); 
			system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Set_$B/Input.txt"));
			
			# Setup QSUB scripts
			$s1_char='\${que}';
			$s2_char='\${output}';
			$s3_char='\${error}';
			system(sprintf("echo \" #PBS -u $UserID\"  >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B"));
			system(sprintf("echo \" #! /bin/bash\"  >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B"));
			system(sprintf("echo \" export HOME=\\\"$OutputDir/workdir$set/\\\"         \"  >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B")); 
			system(sprintf("echo \" export X509_USER_PROXY=\\\"$Proxy\\\"    \"  >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B"));
			system(sprintf("echo \" que=\\\"$QsubQue\\\"\" >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B")); 
			system(sprintf("echo \" output=\\\"Set_$B.qsub.o  \\\"\" >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B")); 
			system(sprintf("echo \" error=\\\"Set_$B.qsub.e  \\\"\" >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B")); 
			system(sprintf("echo \" qsub -q  $s1_char  -o $s2_char -e $s3_char Set_$B.sh\" >> $OutputDir/workdir$set/Set_$B/Qsub_Set_$B"));		
		    }
		    ($a,$b,$c)=split('/',$file);
		    $myfile=$file;
		    if($a =~ m/root/){
			$myfile=$a;
		    }
		    if($b =~ m/root/){
                        $myfile=$b;
			system(sprintf("echo \"notification = Error        \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
                    }
		    if($c =~ m/root/){
                        $myfile=$c;
                    }
		    $myfiletrunc = $myfile;
		    my @wholepath = split /\//, $file;
		    foreach $TreeName (@wholepath){
			if($TreeName =~ m/root/){
			    $myfiletrunc=$TreeName
			}
		    }
#		    system(sprintf("echo \"xrdcp root://sbgse1.in2p3.fr/$file . \"  >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
		    #system(sprintf("echo \"rfcp /dpm/in2p3.fr/home/cms/phedex/$file . \"  >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));  ##rfcp /dpm/in2p3.fr/home/cms/phedex/store/user/cherepan
		    
		    #system(sprintf("echo \"File:  gfal:srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/${DS}0000/$myfiletrunc \"     >> $OutputDir/workdir$set/Set_$B/Input.txt")) ;
			system(sprintf("echo \"File:  root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/$file \"     >> $OutputDir/workdir$set/Set_$B/Input.txt")) ;
		    
                    #system(sprintf("echo \"rm -rf $RemoteDir/$myfiletrunc \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));
		    $A++;
		}
	    }
	}
    }

    printf("Setting up root....\n");
    if($buildRoot==1){
        printf("Building custom root $buildrootversion... Enjoy your coffee!!! ");
	system(sprintf("wget ftp://root.cern.ch/root/root_v$buildrootversion.source.tar.gz"));
	system(sprintf("gzip -dc root_v$buildrootversion.source.tar.gz | tar -xf -"));
	system(sprintf("mkdir $OutputDir/workdir$set/root"));
	system(sprintf("cd root_v$buildrootversion; ./configure --enable-python --enable-roofit --enable-minuit2 --disable-xrootd --disable-sqlite --disable-python --disable-mysql --prefix=$OutputDir/workdir$set/root; make & make install "));
    }

    # Finish Submit script
    system(sprintf("echo \"cd  $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit"));

    # print Instructions
    printf("\n\nInstructions");
    printf("\nPlease make sure you have run:");
#    printf("\nvoms-proxy-init");
    printf("\ngrid-proxy-init"); 
 #   printf("\ngit config --global credential.helper 'cache --timeout=3600'");
 #   printf("\ngit config --global credential.helper cache");
    printf("\nNow you can run the analysis using dcache.");
    printf("\nTo go to the Test workdir: cd  $OutputDir/workdir$set ");
    printf("\nTo compile the code in the workdir: source compile  $UserDir $tauspinner $svfit");
#    printf("\nTo compile the code in the workdir: source compile --useRoot $OutputDir/workdir$set/root/ $UserDir $tauspinner $svfit");
    printf("\nTo submit jobs to the batch queue: source Submit ");
    printf("\nTo combine jobs submitted to the batch queue: source Combine \n");
#    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile  --useRoot $OutputDir/workdir$set/root/ $UserDir; cd $OutputDir/workdir$set/Set_1; source Set_1 | tee log; cd ..\n");
    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile   $UserDir $tauspinner $svfit; cd $OutputDir/workdir$set/Set_1; source Set_1 | tee log; cd ..\n");
} 

if( $ARGV[0] eq "--GRID" ){

    if( $hasroot==0 ){
	printf("ERROR no root version defined ABORTING!!! Please read: ./todo.pl --help \n");
	exit(0);
    }

    $TempDataSetFile=$ARGV[2];
    # Print out input parameters
    printf("Active directory will be: $OutputDir/workdir$set \n");
    printf("Code Repository is:       $CodeDir \n");

    # Open ListofFile.txt
    @DataSets;
    open(DAT, $TempDataSetFile) || die("Could not open file $TempDataSetFile!");
    while ($item = <DAT>) {
	chomp($item);
	print("File: $item \n");
	push(@DataSets,$item);
    }
    close(DAT);

    # Clean Directory in case workdir$set exists
    printf("Cleaning Directories \n");
    system(sprintf("cd $OutputDir"));
    system(sprintf("if [ -d  $OutputDir/workdir$set/ ]; then \n rm -rf $OutputDir/workdir$set \n fi "));
    system(sprintf("mkdir $OutputDir/workdir$set ")); 
    printf("Cleaning complete \n");
    
    #create directory stucture
    system(sprintf("cp $dir/subs $OutputDir/workdir$set "));
    system(sprintf("mkdir  $OutputDir/workdir$set/Code "));
    system(sprintf("mkdir  $OutputDir/workdir$set/Code/i386_linux "));
    system(sprintf("cp -r $CodeDir/* $OutputDir/workdir$set/Code/ "));
    system(sprintf("mkdir $OutputDir/workdir$set/EPS "));
    system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/InputData "));

    # Setup file to retrieve jobs
    system(sprintf("cp $dir/CheckandGet.sh $OutputDir/workdir$set/"));
    system(sprintf("cd $OutputDir/workdir$set/; $dir/subs USERNAME $UserIDCern CheckandGet.sh; cd $dir"));
    system(sprintf("cd $OutputDir/workdir$set/; $dir/subs WORKDIR workdir$set CheckandGet.sh; cd $dir"));
    system(sprintf("cp $dir/Purge_Jobs.sh $OutputDir/workdir$set/"));
    system(sprintf("cp $dir/Cancel_Jobs.sh $OutputDir/workdir$set/"));
    system(sprintf("cp $dir/Run.sh $OutputDir/workdir$set/"));
    system(sprintf("cp $dir/ResubmitFailedJobs.py $OutputDir/workdir$set/"));

    # generate compile script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"source config \\\$\@ \"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake all \" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/compile")) ;
 
    # Generate Combine script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Combine")) ;
    system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Combine"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Combine")) ; 
    system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe \" >> $OutputDir/workdir$set/Combine")) ;

    # Generate Combine Input
    system(sprintf("cp $InputFile $OutputDir/workdir$set/Input.txt "));
    system(sprintf("cd $OutputDir/workdir$set; subs '{SET}' COMBINE $OutputDir/workdir$set/Input.txt; cd $dir "));
    system(sprintf("cd $OutputDir/workdir$set; subs '{FileDir}' COMBINE $OutputDir/workdir$set/Input.txt; cd $dir "));
    system(sprintf("echo \"Mode: RECONSTRUCT\" >> $OutputDir/workdir$set/Input.txt"));
    system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Input.txt"));

    # Setup Condor Combine scripts
    system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"executable   = Combine      \"  >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"output       = Combine-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"error        = Combine-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"log          = Combine-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Condor_Combine")); 			
    system(sprintf("echo \"queue = 1 \" >> $OutputDir/workdir$set/Condor_Combine"));

    # Start Submit script
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Submit")) ; 
# This script run the grid jobs and then combines the output
    system(sprintf("echo 'if [ \"\${1}\"  == \"--help\" ] ; then ' >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo '  echo \"Script to submit grid jobs. \" ' >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo '  echo \"Options for running this this script\" ' >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo '  echo \"source Submit --Setup           Creates new tarball and ships it to the grid. Cleans history of submitted jobs. \" ' >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo '  echo \"source Submit --Submit          Submits all unsubmitted jobs using installed tarball. \" ' >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo '  echo \"source Submit --SetupAndSubmit  Runs --Setup and then --Submit. \" ' >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo 'fi ' >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"verbosity=\\\$(grep SetLevel Code/Analysis.cxx | grep -c -e Debug -e Verbose)\" >> $OutputDir/workdir$set/Submit")) ;
	system(sprintf("echo \"  if [[ \\\${verbosity} -ne 0 ]]; then \" >> $OutputDir/workdir$set/Submit")) ;
	system(sprintf("echo \"    echo 'ERROR: Please make sure to set the verbosity level to Info in Analysis.cxx, otherwise your log-files will break DCache! Abort...' \" >> $OutputDir/workdir$set/Submit"));
	system(sprintf("echo \"    return \\\${verbosity}\" >> $OutputDir/workdir$set/Submit"));
	system(sprintf("echo \"  fi\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo 'if [ \"\${1}\" == \"--Setup\" ] || [ \"\${1}\" == \"--SetupAndSubmit\" ]; then ' >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"  cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"  if [ -f $OutputDir/workdir$set/Set*/out ]; then \n    rm $OutputDir/workdir$set/Set*/out;\n    rm $OutputDir/workdir$set/Set*/err;\n    rm $OutputDir/workdir$set/Set*/*.tar; \n  fi  \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  echo 'Creating tarballs and installing on the GRID... '\" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"  if [ -f workdir$set.tar  ]; then \n    rm workdir$set.tar \n  fi \n  tar -cf workdir$set.tar root Code;\" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"  isSkim=\\\$(srmls -recursion_depth=2 srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/ | grep workdir$set.tar  | wc -l) \n  if [ \\\$isSkim == 1 ]; then \" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"    srmrm  srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set.tar \n  fi \" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"  srmcp  file:////$OutputDir/workdir$set/workdir$set.tar srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set.tar \" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"  echo 'Complete '\" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"  if [ -f $OutputDir/workdir$set/jobs_submitted ]; then \n    rm $OutputDir/workdir$set/jobs_submitted \n  fi \n  touch $OutputDir/workdir$set/jobs_submitted\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  if [ -f $OutputDir/workdir$set/jobs_complete ]; then  \n  rm $OutputDir/workdir$set/jobs_complete  \n  fi \n  touch $OutputDir/workdir$set/jobs_submitted\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  if [ -f $OutputDir/workdir$set/jobs_submittedOrComplete ]; then  \n    rm $OutputDir/workdir$set/jobs_submittedOrComplete  \n  fi \n  touch $OutputDir/workdir$set/jobs_submittedOrComplete\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  touch jobs_log \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  echo 'Configuring storage on the GRID...' \" >> $OutputDir/workdir$set/Submit"));
 system(sprintf("echo \"  hasdir=\\\$(srmls -recursion_depth=2 srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/ | grep /pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/  | wc -l) \n  if [ \\\$hasdir == 0 ]; then \" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"    srmmkdir srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set \n  fi \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo 'fi ' >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("chmod u+x $OutputDir/workdir$set/Submit")) ;

    # get file list from dcache to remove old skims
    system(sprintf("echo '\n\n\nif [ \"\${1}\" == \"--Submit\" ] || [ \"\${1}\" == \"--SetupAndSubmit\" ]; then ' >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"  cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"  #find list of Skim file from previous job\n  touch jobs_submittedOrComplete \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"  touch listOfSrmlsFileNames \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"  filecount=0 \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"  readfiles=1 \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"  while [ \\\$readfiles -gt \\\$filecount ]; do \" >> $OutputDir/workdir$set/Submit")) ;
     system(sprintf("echo \"    srmls -count=999 -offset=\\\$filecount -recursion_depth=2 srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/ >& listOfSrmlsFileNames \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"    let filecount=filecount+999 \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"    readfiles=\\`cat listOfSrmlsFileNames | wc -l\\`; \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"  done \" >> $OutputDir/workdir$set/Submit")) ;
    
    


    $B=0;
    for($l=0;$l<2; $l++){
	system(sprintf("echo \"notification = Error \" >> $OutputDir/workdir$set/Condor_Combine"));
	$max=1;
	foreach $DS (@DataSets){
	    if(($l==0 && ($DS =~ m/Data/)) || ($l==1 && !($DS =~ m/Data/))){
		if($l==0){
		    $max=$maxdata;
		}
		else{
		    $max=$maxmc;
			if($DS =~ m/embed/){
			    $max=$maxemb;
			}
		}
		printf("\n\nStarting Loop $l \n");
		$A=$maxdata+$maxmc+$maxemb+10;

		# find the root files for the current DataSet (DS)
		printf("Accessing Directory  $DS \n");
		system(sprintf("touch junk"));

		system(sprintf("touch junk"));
                system(sprintf("uberftp $ubergridsite \"cd $DS; ls */  \" | grep root >& junk "));
                system(sprintf("cat junk | awk '{print \$8}' >& junk1"));

		# Get list of files in dcache dir
		@files=(); 
		open(DAT, "junk1");
		while ($item = <DAT>) {
		    chomp($item);
		    push(@files,$item);
		}
		close(DAT);

		system(sprintf("rm junk"));
		system(sprintf("rm junk1"));
		$nfiles = @files;
		$idx=0;
		
		
		
		foreach $file (@files){
		    $idx++;
		    #printf("$DS/$file Set= $B  Index=  $A   Max.= $max N Files = $nfiles Current File = $idx \n");
		    if($A > $max ){
			$A=1;
			$B++;
			
                        # Add Set information to Combining scripts and Input.txt
			system(sprintf("echo \"File: $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Input.txt ")) ;

			# Create and configure Set_$B dir
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B ")) ;
			system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/Set_$B/InputData "));
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B/EPS ")) ;
			
			# Setup Set_$B.sh
			system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("echo \"echo 'Starting Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"export workdir=\\\"$OutputDir/workdir$set/\\\"\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"source $OutputDir/workdir$set/Set_$B/Set_$B-get.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ; 
			system(sprintf("chmod +x $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"mkdir  /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cp -r *   /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe 2>&1 | tee >(sed -r \\\"s/\\\\x1B\\\\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g\\\" > Set_$B.output) \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cp -r *  $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			#system(sprintf("echo \"source $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"rm -r  /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));			
			system(sprintf("echo \"echo 'Completed Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;

            # Setup Set_$B-GRID.sh 
            system(sprintf("echo \"echo 'Starting Job'; ls; echo \\\$PWD \" > $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("echo 'srmcp  srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set.tar file:////\$PWD/workdir$set.tar ; tar -xf workdir$set.tar; ' >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh  "));
			system(sprintf("echo \"cp \\\$PWD/Inputgrid.txt \\\$PWD/Input.txt\" >>  $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh  "));
			system(sprintf("echo \"cp -r \\\$PWD/Code/InputData/ \\\$PWD/InputData \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("chmod +x $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			#system(sprintf("echo \"source \\\$PWD/Set_$B-getGRID.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
            system(sprintf("echo \"cd Code/; chmod +x  \\\$PWD/config  \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
            system(sprintf("echo \"echo \\\$PWD | awk '{split(\\\$1,a,\\\"Code\\\"); print \\\"#! /bin/bash  \\n source \\\" \\\$1 \\\"/config --useRoot \\\" a[1] \\\"root/ \\\"}'  > \\\$PWD/junk; source \\\$PWD/junk; \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("echo \" cd ../; echo \\\$PWD; ls  \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			#system(sprintf("echo \"chmod +x \\\$PWD/subs; \\\$PWD/subs '/user/scratch/$UserID' \\\$PWD Input.txt \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("echo \"export LD_LIBRARY_PATH=\\\$LD_LIBRARY_PATH:\\\$ROOTSYS/lib/\" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("echo \"export LHAPATH=\\\$PWD/Code/TauSpiner/lhapdf/share/lhapdf/PDFsets/\" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("echo \"echo 'System Configured.'; printenv; ls ; echo \\\$PWD ; \\\$PWD/Code/Analysis.exe | sed -r \\\"s/\\\\x1B\\\\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g\\\" \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			#
			system(sprintf("echo \"if [ -f SKIMMED_NTUP.root ]; then \n \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh "));
			system(sprintf("echo \"srmcp  file:////\\\$PWD/SKIMMED_NTUP.root srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/SKIMMED_NTUP_$B.root \n fi \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh "));
			#
			system(sprintf("echo \"echo 'Completed Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
                        #
			system(sprintf("echo \"#! /bin/bash \" >> $OutputDir/workdir$set/Set_$B/GRIDRetrieve.sh "));
			system(sprintf("echo \"if [ -f $OutputDir/workdir$set/fileListOnGrid ]; then \" >> $OutputDir/workdir$set/Set_$B/GRIDRetrieve.sh "));
			system(sprintf("echo \"isSkim=\\\$(cat $OutputDir/workdir$set/fileListOnGrid | grep SKIMMED_NTUP_$B.root  | wc -l) \" >> $OutputDir/workdir$set/Set_$B/GRIDRetrieve.sh "));
			system(sprintf("echo \"else \" >> $OutputDir/workdir$set/Set_$B/GRIDRetrieve.sh "));
    		system(sprintf("echo \"isSkim=\\\$(srmls -recursion_depth=2 srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/ | grep SKIMMED_NTUP_$B.root  | wc -l) \" >> $OutputDir/workdir$set/Set_$B/GRIDRetrieve.sh "));
			system(sprintf("echo \"fi \" >> $OutputDir/workdir$set/Set_$B/GRIDRetrieve.sh "));
			system(sprintf("echo \"if [ \\\$isSkim == 1 ]; then \" >> $OutputDir/workdir$set/Set_$B/GRIDRetrieve.sh "));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B; srmcp  srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/SKIMMED_NTUP_$B.root  file:////\\\$PWD/SKIMMED_NTUP.root; cd $OutputDir/workdir$set/ \n fi \" >> $OutputDir/workdir$set/Set_$B/GRIDRetrieve.sh "));
			
                        # Setup Set_$B_get.sh Set_$B-getGRID and Set_$B_clean.sh
			#system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			#system(sprintf("echo \"mkdir /user/scratch/$UserID \" >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			#system(sprintf("echo \"cd /user/scratch/$UserID \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh")); 
                        #system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));
                        #system(sprintf("echo \"cd /user/scratch/$UserID \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));
                        #system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-getGRID.sh"));

			# Setup Input.txt
			system(sprintf("cp   $InputFile $OutputDir/workdir$set/Set_$B/Input.txt ")); 
			system(sprintf("cd  $OutputDir/workdir$set/Set_$B; $dir/subs '{SET}' Set_$B Input.txt; cd $dir "));
			system(sprintf("cd  $OutputDir/workdir$set/Set_$B; $dir/subs '{FileDir}' $DS Input.txt; cd $dir "));
			system(sprintf("echo \"Mode: ANALYSIS\" >> $OutputDir/workdir$set/Set_$B/Input.txt")); 
			system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Set_$B/Input.txt"));
			system(sprintf("cp   $InputFile $OutputDir/workdir$set/Set_$B/Inputgrid.txt "));
			system(sprintf("cd  $OutputDir/workdir$set/Set_$B; $dir/subs '{SET}' Set_$B Inputgrid.txt; cd $dir "));
			system(sprintf("cd  $OutputDir/workdir$set/Set_$B; $dir/subs '{FileDir}' $DS Inputgrid.txt; cd $dir "));
                        system(sprintf("echo \"Mode: ANALYSIS\" >> $OutputDir/workdir$set/Set_$B/Inputgrid.txt"));
                        system(sprintf("echo \"RunType: GRID\" >> $OutputDir/workdir$set/Set_$B/Inputgrid.txt"));

			# Setup .jdl file
                        system(sprintf("echo '['  >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo 'JobType = \"Normal\";'  >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
                        system(sprintf("echo 'Executable    = \"Set_$B-GRID.sh\"; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo 'Arguments     = \" \"; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo 'StdOutput     = \"out\"; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo 'StdError      = \"err\"; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo 'InputSandbox  = {\"Set_$B/Set_$B-GRID.sh\",\"subs\",\"Set_$B/Inputgrid.txt\"}; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo 'OutputSandbox = {\"out\",\"err\",\"*_ANALYSIS_*.root\"}; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo 'Requirements  = (RegExp(\"rwth-aachen.de\", other.GlueCEUniqueId)) && (RegExp(\"cream\", other.GlueCEUniqueId)); ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo 'OutputSandboxBaseDestUri=\"gsiftp://localhost\" '  >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			#system(sprintf("echo ' Retry         = 0; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl")); 
			system(sprintf("echo ']' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl")); 
			
		    }
		    ($x,$y)=split($DS,$file);
		    if($x =~ m/root/){
			$y=$x;
		    }
		    ($a,$b,$c,$d,$e)=split('/',$y);
		    $file=$y;
		    if($a =~ m/root/){
			$myfile=$a;
		    }
		    if($b =~ m/root/){
                        $myfile=$b;
                    }
		    if($c =~ m/root/){
                        $myfile=$c;
                    }
                    if($d=~ m/root/){
                        $myfile=$d;
                    }

                    if($e =~ m/root/){
                        $myfile=$e;
                    }

		    system(sprintf("echo \"dccp  dcap://grid-dcap-extern.physik.rwth-aachen.de/$DS/$file . \"  >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
		    system(sprintf("echo \"File:  /user/scratch/$UserID/$myfile \"     >> $OutputDir/workdir$set/Set_$B/Input.txt")) ;
		    system(sprintf("echo \"File:  dcap://$dcapgridsite/$DS/$myfile \"     >> $OutputDir/workdir$set/Set_$B/Inputgrid.txt")) ;
		    $A++;
		}
	    }
	}
    }


    printf("Setting up root....\n");
    if($buildRoot==1){
        printf("Building custom root $buildrootversion... Enjoy your coffee!!! ");
	system(sprintf("wget ftp://root.cern.ch/root/root_v$buildrootversion.source.tar.gz"));
	system(sprintf("gzip -dc root_v$buildrootversion.source.tar.gz | tar -xf -"));
	system(sprintf("mkdir $OutputDir/workdir$set/root"));
	system(sprintf("cd root_v$buildrootversion; ./configure --enable-python --enable-roofit --enable-minuit2 --disable-xrootd --disable-sqlite --disable-python  --disable-mysql --prefix=$OutputDir/workdir$set/root; make & make install "));
    }
    else{
        printf("Copying local root $MYROOTSYS ");
        system(sprintf("mkdir $OutputDir/workdir$set/root/"));
        system(sprintf("cp -r $MYROOTSYS/* $OutputDir/workdir$set/root/"));
    }
    
    # Finish Submit script
    system(sprintf("echo \"\n\n# Now setup up submission (try 3 times for good luck) \n  for (( l=1; l<=3; l++ )); do \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"    for (( c=1; c<=$B; c++ )); do \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"      issub=\\\$(grep \\\"workdir$set/Set_\\\${c} \\\" $OutputDir/workdir$set/jobs_submittedOrComplete | wc -l) \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"      if [[ \\\${issub} -eq 0 ]]; then  \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"        cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"        echo 'Submitting Set_'\\\${c}'/GRIDJob.jdl' \" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"        isSkim=\\\$(cat listOfSrmlsFileNames | grep SKIMMED_NTUP_\\\${c}.root | wc -l) \n      if [[ \\\${isSkim} -eq 1 ]]; then \" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"          srmrm srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/SKIMMED_NTUP_\\\${c}.root \n     fi \" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"          glite-ce-job-submit -a -r grid-ce.physik.rwth-aachen.de:8443/$Queue $OutputDir/workdir$set/Set_\\\${c}/GRIDJob.jdl | tee junk ; cat junk >> jobs_log; cat junk | grep https | awk -v idx=\\\${c} '{print \\\$1 \\\" $OutputDir/workdir$set/Set_\\\" idx \\\" \\\"}' | tee -a $OutputDir/workdir$set/jobs_submittedOrComplete >> $OutputDir/workdir$set/jobs_submitted ; rm junk \" >> $OutputDir/workdir$set/Submit"));
    #system(sprintf("echo \"          cp $OutputDir/workdir$set/jobs_submitted $OutputDir/workdir$set/jobs_submittedOrComplete \" >> $OutputDir/workdir$set/Submit")) ; 
    system(sprintf("echo \"      fi \n    done\n  done \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  njobs=\\\$(cat jobs_submittedOrComplete | wc -l)  \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  if [[ \\\${njobs} -ne $B ]]; then \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"    echo \\\${njobs} ' out of $B were submitted. Check that the grid is not experiencing technical difficulties and retry: source Submit --Submit ' \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  fi\" >> $OutputDir/workdir$set/Submit"));

    #Clean up and Summary
    system(sprintf("echo \"\n\n\n#Summary and Cleanup\n  cd  $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  rm  listOfSrmlsFileNames \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  echo 'number of jobs submitted:'\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  cat jobs_submitted | grep -c Set_\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  echo 'total number of sets in this directory:'\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"  ls -l | grep -c Set_\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"fi\" >> $OutputDir/workdir$set/Submit"));

    # print Instructions
    printf("\n\nInstructions");
    printf("\nPlease make sure you have run:");
    printf("\nvoms-proxy-init");
    printf("\ngrid-proxy-init");
    printf("\ngit config --global credential.helper 'cache --timeout=3600'");
    printf("\ngit config --global credential.helper cache");
    printf("\nNow you can run the analysis using dcache.");
    printf("\nTo go to the Test workdir: cd  $OutputDir/workdir$set ");
#    printf("\nTo compile the code in the workdir: source compile --useRoot $OutputDir/workdir$set/root/ $UserDir $tauspinner $svfit");
    printf("\nTo compile the code in the workdir: source compile  $UserDir $tauspinner $svfit");
    printf("\nTo submit jobs to the GRID: source Submit ");
    printf("\nTo check the status of the GRID jobs and download finished jobs: source CheckandGet.sh");
    printf("\nTo additionally print the details of all the jobs: source CheckandGet.sh  --detailed");
#    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile  --useRoot $OutputDir/workdir$set/root/ $UserDir; cd $OutputDir/workdir$set/Set_1; source Set_1.sh; cd ..\n");
    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile    $UserDir $tauspinner $svfit; cd $OutputDir/workdir$set/Set_1; source Set_1.sh; cd ..\n");
    
} 
