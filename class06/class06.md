# Class06
Anyoleth Alarcon (A17347293)

- [1. Function Basics](#1-function-basics)
- [2. Generate DNA Sequence](#2-generate-dna-sequence)
- [3. Generate Protein Function](#3-generate-protein-function)

## 1. Function Basics

Let’s start writing out our first silly function to add some numbers:

Every R function has 3 components:

- name (we get to pick this)
- input arguments (there can be loads of these separated by comma)
- the body (the R code that does the work)

``` r
add <- function(x, y=10, z=0){
  x + y + z
}
```

I can just use this function like any other function as long as R knows
about it (i.e. run the code chunk)

``` r
add(1, 100)
```

    [1] 101

``` r
add(x=c(1,2,3,4), y=100)
```

    [1] 101 102 103 104

``` r
add(1)
```

    [1] 11

Functions can have “required” input arguments and “optional” input
arguments. The optional arguments are defined with an equals default
value (`y=10`) in the function definition.

``` r
add(1,100,10)
```

    [1] 111

> Q. Write a function to return a DNA sequence of a user specified
> length? Call it `generate_DNA()`

The `sample()` function can help here

``` r
#generate_DNA <- function(size=5) {}

students <- c("jeff", "jeremy", "peter")

sample(students, size=5, replace=TRUE)
```

    [1] "peter"  "jeremy" "jeremy" "jeremy" "peter" 

Now work with `bases` rather than `students`

``` r
bases <- c("A", "C", "G", "T")
sample(bases, size=10, replace=TRUE)
```

     [1] "A" "A" "A" "C" "G" "A" "T" "C" "A" "A"

## 2. Generate DNA Sequence

Now I have a working `snippet` of code I can use this as the body of my
first function version here:

``` r
generate_DNA <- function(size=5) {
  bases <- c("A", "C", "G", "T")
  sample(bases, size=size, replace=TRUE)
}
```

``` r
generate_DNA(100)
```

      [1] "C" "T" "T" "C" "G" "T" "A" "C" "T" "G" "C" "G" "T" "A" "A" "G" "G" "T"
     [19] "G" "G" "T" "T" "T" "A" "G" "G" "C" "T" "T" "T" "T" "T" "G" "G" "T" "C"
     [37] "T" "T" "A" "T" "G" "G" "C" "T" "A" "G" "T" "C" "G" "T" "T" "C" "G" "C"
     [55] "A" "A" "A" "T" "T" "T" "T" "C" "G" "C" "T" "G" "G" "C" "T" "G" "G" "A"
     [73] "C" "A" "T" "G" "C" "C" "G" "C" "G" "C" "C" "C" "C" "A" "T" "T" "C" "G"
     [91] "T" "T" "C" "G" "A" "C" "A" "G" "G" "A"

``` r
generate_DNA(100)
```

      [1] "G" "A" "A" "T" "G" "T" "C" "C" "G" "G" "C" "C" "G" "T" "C" "C" "A" "C"
     [19] "C" "C" "C" "T" "C" "C" "T" "A" "C" "T" "T" "A" "G" "A" "A" "T" "T" "G"
     [37] "T" "C" "T" "A" "G" "C" "C" "G" "C" "T" "G" "G" "C" "G" "T" "A" "G" "T"
     [55] "T" "T" "C" "T" "T" "G" "T" "T" "T" "G" "C" "T" "G" "T" "G" "A" "C" "A"
     [73] "C" "T" "C" "C" "C" "A" "G" "T" "G" "T" "G" "C" "A" "C" "A" "A" "G" "G"
     [91] "C" "A" "T" "C" "C" "T" "G" "C" "A" "C"

I want the ability to return a sequence like “AGTACCTG” i.e. a one
element vector where the bases are all together.

``` r
generate_DNA <- function(size=5, together=TRUE) {
  bases <- c("A", "C", "G", "T")
  sequence <- sample(bases, size=size, replace=TRUE)
  if(together){
  sequence <- paste(sequence, collapse = "")
    
  }
  return(sequence)
}
```

``` r
generate_DNA(together=FALSE)
```

    [1] "C" "G" "G" "C" "A"

``` r
generate_DNA(20)
```

    [1] "CGTGTAGGAAGTTGATTTAA"

## 3. Generate Protein Function

> Q. Write a protein sequence generating function that will return
> sequences of a user specified length?

We can get the set of 20 natural amino-acids from the **bio3d** package.

``` r
aa <- bio3d::aa.table$aa1[1:20]
```

``` r
generate_protein <- function(size=6, together=TRUE){
  
  ##Get the 20 amino-acids as a vector
  aa <- bio3d::aa.table$aa1[1:20]
  sequence <- sample(aa, size=size, replace=TRUE)
  
  ##Optionally return a single element string
  if(together){
  sequence <- paste(sequence, collapse = "")
    
  }
  return(sequence)
}
```

``` r
generate_protein(together = F)
```

    [1] "N" "M" "G" "Y" "I" "G"

``` r
generate_protein(7)
```

    [1] "IAWGPWD"

> Q. Generate random protein sequences of length 6 to 12 amino acids.

``` r
##Generate_protein(size=6:12)
```

We can fix this inability to generate multiple sequences by either
editing and adding the function body code (e.g. a for loop) or by using
the R **apply** family of utility functions.

``` r
sapply(6:12, generate_protein)
```

    [1] "SLNYLV"       "SNLVEIV"      "YCNIYHWD"     "LQVLGNHEG"    "CMFFCWMNDR"  
    [6] "QFIWNQVCNWC"  "YDKWMDGNTCDE"

It would be cool and useful if I could get FASTA format output

``` r
ans <- sapply(6:12, generate_protein)
ans
```

    [1] "RMGGKP"       "YRGAFIK"      "RFPICRWM"     "MTWGVWDSW"    "DKKDENFEAY"  
    [6] "NWAFADWTFKD"  "LIPHIGTSWCKS"

``` r
cat(ans, sep="\n")
```

    RMGGKP
    YRGAFIK
    RFPICRWM
    MTWGVWDSW
    DKKDENFEAY
    NWAFADWTFKD
    LIPHIGTSWCKS

I want this to look like

    >ID.6
    PADREN
    >ID.7
    RTNVGPT
    >ID.8
    YNGFMNYF
    MECRCGCCW
    VCDDMGTDHN
    WQDKVHTAYDA
    WSHAQVWWKGLT
    and so forth

The functions `paste()` and `cat()` can help us here…

``` r
paste(">ID.", 7:12, sep = "")
```

    [1] ">ID.7"  ">ID.8"  ">ID.9"  ">ID.10" ">ID.11" ">ID.12"

``` r
cat( paste(">ID.", 6:12,"\n", ans, sep = ""), sep="\n")
```

    >ID.6
    RMGGKP
    >ID.7
    YRGAFIK
    >ID.8
    RFPICRWM
    >ID.9
    MTWGVWDSW
    >ID.10
    DKKDENFEAY
    >ID.11
    NWAFADWTFKD
    >ID.12
    LIPHIGTSWCKS

``` r
id.line <- paste(">ID.",6:12, sep="")
id.line
```

    [1] ">ID.6"  ">ID.7"  ">ID.8"  ">ID.9"  ">ID.10" ">ID.11" ">ID.12"

``` r
seq.line <- paste(id.line, ans, sep="\n")
cat(seq.line, sep="\n")
```

    >ID.6
    RMGGKP
    >ID.7
    YRGAFIK
    >ID.8
    RFPICRWM
    >ID.9
    MTWGVWDSW
    >ID.10
    DKKDENFEAY
    >ID.11
    NWAFADWTFKD
    >ID.12
    LIPHIGTSWCKS

> Q. Determine if these sequences can be found in nature or are they
> unique?

I BLASTp searched my FASTA format sequences against NR and found that
length 6 and 7 are not unique and can be found in the databases with
100% coverage and 100% identity.

Random sequence of length of 8 and above are unique and cannot be found
in the databases.
