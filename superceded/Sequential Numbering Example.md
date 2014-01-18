# Incremental Figure Numbers


```r
fn = local({
    i = 0
    function(x) {
        i <<- i + 1
        paste("Figure ", i, ": ", x, sep = "")
    }
})
```



```r
plot(1:10)
```

![Figure 1: this is one plot](figure/test-a.png) 



```r
plot(rnorm(10))
```

![Figure 2: another plot](figure/test-b.png) 

