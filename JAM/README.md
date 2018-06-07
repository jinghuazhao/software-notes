## Setup

The package is available from https://github.com/pjnewcombe/R2BGLiMS.

Note that JAM requires Java 1.8 so call to Java -jar inside the function needs to
reflect this, not straightforward with `install_github()` from `devtools` but one needs to
clone the package, modify the R source code and then install,
```
git clone https://github.com/pjnewcombe/R2BGLiMS
### change java to java-1.8 in R2BGLiMS/R/R2BGLiMS.R
R CMD INSTALL R2BGLiMS
```

## Compiling

The information is unavailable from the documentation, but at least can be achieved this with [netbeans](https://netbeans.org/).

## Old notes
```bash
21-7-2017 MRC-Epid JHZ

export BGLiMS=/genetics/bin/BGLiMS
export JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.131-0.b11.el6_9.x86_64

# Jama

cd $BGLiMS
$JAVA_HOME/bin/javac Jama/util/*java
$JAVA_HOME/bin/javac -classpath $BGLiMS Jama/*java
$JAVA_HOME/bin/jar cvf Jama.jar Jama/*class Jama/util/*class

# ssj

export SSJHOME=ssj
export LD_LIBRARY_PATH=$SSJHOME/lib:$LD_LIBRARY_PATH
export CLASSPATH=.:$SSJHOME/lib/ssj.jar:$SSJHOME/lib/colt.jar:$SSJHOME/lib/tcode.jar:$CLASSPATH

# MyClassLibrary

cd $BGLiMS/MyClassLibrary/src
ln -sf $BGLiMS/Jama
$JAVA_HOME/bin/javac Methods/*java
$JAVA_HOME/bin/javac Objects/*java
$JAVA_HOME/bin/jar cvf MyClassLibrary.jar Methods/*class Objects/*class

# BGLiMS.jar

cd $BGLiMS/src
$JAVA_HOME/bin/javac bglims/*java
$JAVA_HOME/bin/jar cvf bglims.jar bglims/*class

cd $BGLiMS
ln -sf $BGLiMS/src/bglims.jar
ln -sf MyClassLibrary/src/MyClassLibrary.jar
# $JAVA_HOME/bin/jar cvf BGLiMS.jar bglims.jar MyClassLibrary.jar
```
