Bootstrap: docker
From: ubuntu:bionic

%runscript
echo "Running the container.."

%post

# Install latest version of R (https://linuxize.com/post/how-to-install-r-on-ubuntu-18-04/)
#  (do it this way otherwise by default will install 3.2)

apt-get -y update

apt-get install -y apt-transport-https software-properties-common
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

apt-get -y update
apt-get -y install r-base-core

apt-get -y install tk-dev

# Install required R packages
R --slave -e 'install.packages("RColorBrewer", dependencies=TRUE, repos = "http://cran.us.r-project.org")'
R --slave -e 'install.packages("TeachingDemos", dependencies=TRUE, repos = "http://cran.us.r-project.org")'




# Installing perl modules

cpan -P File::chdir
cpan -P Config::General
cpan -P Data::Dumper
cpan -P Excel::Writer::XLSX
cpan -P DBI
cpan -P DBD::SQLite



# Installing LaTeX stuff
#apt-get -y --no-install-recommends install texlive-latex-base texlive-fonts-recommended texlive-latex-extra lmodern
apt-get -y install texlive-latex-base texlive-fonts-recommended texlive-latex-extra lmodern
apt-get -y install wget


# Installing extra software
apt-get -y install libssl-dev
apt-get -y install libcurl4-openssl-dev
apt-get -y install git


# Getting GenomeChronicler Repo
git clone git@github.com:PGP-UK/GenomeChronicler.git


# Running GenomeChronicler Setup
cd GenomeChronicler
bash SetupMeFirst.sh
ln -sf GenomeChronicler_mainDruid.pl genomechronicler
chmod +x GenomeChronicler_mainDruid.pl
chmod +x genomechronicler


%environment
export LC_ALL=C
export PATH=$PATH:$PWD

#genomechronicler


