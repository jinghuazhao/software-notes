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
