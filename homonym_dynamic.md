Homonyms, French Names, dynamics
================
Arthur Charpentier & Baptiste Coulmont
July 2017

``` r
load("/home/arthur/Dropbox/Dropbox (Télécom SudParis)/Nouveau dossier 4/charpentier.Rdata")
head(patronymes)
```

    ##            nom np1640 np4165 np6690
    ## 1      -IZEROS      0      0      1
    ## 2            A      0      0      6
    ## 3     A ABDOUL      0      0      2
    ## 4 A AHAD HADAD      0      0      2
    ## 5  A AL AZZAWI      0      0      1
    ## 6   A ALFAISAL      0      0      1

``` r
head(prenoms)
```

    ##         preusuel np1640 np4165 np6690
    ## 1 _PRENOMS_RARES 106057 118352 306804
    ## 2              A     NA     NA      3
    ## 3          AADEL     NA     NA     28
    ## 4          AADIL     NA     NA    113
    ## 5          AARON      6     22    146
    ## 6          AAZIZ     NA     NA     17

``` r
dim(patronymes)
```

    ## [1] 1329273       4

``` r
x1=prenoms[,2]
x2=prenoms[,3]
x3=prenoms[,4]
statistiques=function(x,which=NA){
x[is.na(x)]=0
p1=sum(sort(x[-1],decreasing=TRUE)[1])/(sum(x))
p2=sum(sort(x[-1],decreasing=TRUE)[1:10])/(sum(x))
p3=sum(sort(x[-1],decreasing=TRUE)[1:100])/(sum(x))
sv=NA
if(which=="prenom"){
n=x[1]
vn=c(round(n/2),round(n/4),length(x[-1]))
sv=sum(vn)}
if(which=="nom"){
sv=length(x[x>0])}
c(sv/1e6,100*p1,100*p2,100*p3)
}
```

Period *1916-1940*, number of first names (in millions)

``` r
s1=statistiques(x1,which="prenom")
s1[1]
```

    ## [1] 0.095708

use here the following rule for *rare* names : if *n* have a rare names, *n*/2 have a unique one, and for the remaining half, we have pairs, i.e. there are 3*n*/4 first names for those *n* people

Most popular first name (in %)

``` r
s1[2]
```

    ## [1] 5.175006

while the top ten represents (in %)

``` r
s1[3]
```

    ## [1] 25.17757

and the top 100

``` r
s1[4]
```

    ## [1] 79.05222

For *1941-1965* and *1966-1990*, we have (numbers in millions for the first column, and then percentages for the second colum)

``` r
statistiques(x2,which="prenom")
```

    ## [1]  0.10493  3.78186 20.49576 72.61395

``` r
statistiques(x3,which="prenom")
```

    ## [1]  0.246269  1.508596 12.590275 56.975909

For the last name, load the datasets

``` r
y1=patronymes[,2]
y2=patronymes[,3]
y3=patronymes[,4]
```

For the periode *1916-1940*, the number of last names (in millions) is

``` r
s1=statistiques(y1,which="nom")
s1[1]
```

    ## [1] 0.638101

the most popular represent (in %)

``` r
s1[2]
```

    ## [1] 0.3616469

of the top list of last names. For the top 10, (in %)

``` r
s1[3]
```

    ## [1] 1.83409

and for the top 100

``` r
s1[4]
```

    ## [1] 8.663211

For periods *1941-1965* and *1966-1990*, we have (number in millions for the first column, and then percentage for the second one)

``` r
statistiques(y2,which="nom")
```

    ## [1] 0.6693460 0.3518078 1.7607637 8.4118325

``` r
statistiques(y3,which="nom")
```

    ## [1] 0.8140410 0.3155469 1.5723140 7.8381420

The proportion of homonyms can be computed using the following code (first and last names are selected independently)

``` r
prophomonymes=function(n,z){
  x=z$x
  x[is.na(x)]=0
  x0=x[-1]
  listex=c(rep(1,round(x[1]/2)),rep(2,round(x[1]/4)),x0[x0>0])
  sx=sample(1:length(listex),size=n,replace=TRUE,prob=listex)
  y=z$y
  listey=y[y>0]
  sy=sample(1:length(listey),size=n,replace=TRUE,prob=listey)
  sz=paste(sx,sy)
  t=table(sz)
  sum(t[t>1])/n}
sim=function(n,ns=50){
  P=matrix(NA,ns,3)
  for(i in 1:ns){
    P[i,1]=prophomonymes(list(x=x1,y=y1),n=n)
    P[i,2]=prophomonymes(list(x=x2,y=y2),n=n)
    P[i,3]=prophomonymes(list(x=x3,y=y3),n=n)
  }
  apply(P,2,mean)
}
```

Hence, for a group of *n* = 100, 000 people, the porportion of people that have an homonym is (in %), for each periods (the oldest on the left),

``` r
sim(n=100000,ns=50)
```

    ## [1] 0.0940052 0.0734224 0.0404904

with 50 simulations in those Monte Carlo simulations

``` r
P=matrix(NA,5,3)
P[1,]=sim(10000,50)
P[2,]=sim(20000,50)
P[3,]=sim(50000,50)
P[4,]=sim(100000,50)
P[5,]=sim(200000,50)
```

To visualize it, use

``` r
proba=data.frame(annees=rep(c("1916-1940","1941-1965","1966-1990"),5),
                 proba=c(P[1,],P[2,],P[3,],P[4,],P[5,])*100,
                 Taille=as.factor(rep(c(10000,20000,50000,100000,200000),each=3)))
library(ggplot2)
ggplot(proba,aes(annees,proba,color=Taille,group=Taille))+geom_point()+geom_line()+
  scale_y_log10()+ylab("Probabilité (%)")+xlab("") + theme(legend.position="top")
```

![](homonym_dynamic_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-16-1.png)
