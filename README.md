# datasciencecoursera

*** EXAMPLE OF A MARKDOWN DOCUMENT ***

![R icon](https://d1q6f0aelx0por.cloudfront.net/product-logos/87571b58-883f-427b-a52b-b2784a035c6a-r-base.png)

***

    hh <- hclust(dist(dataMatrix))
    dataMatrixOrdered <- dataMatrix[hh$order, ]
    par(mfrow = c(1,3))
    image(1:10,1:40, t(dataMatrixOrdered))[,nrow(dataMatrixOrdered):1] 
    plot(rowMeans(dataMatrixOrdered), 40:1,xlab = "Row Mean", ylab = "Row", pch = 19)
    plot(colMeans(dataMatrixOrdered), xlab = "Column", ylab = "Column mean", pch = 19)

##Link to Markdown tutorial
<https://blog.wax-o.com/2014/04/tutorial-short-guide-to-start-with-markdown/>
