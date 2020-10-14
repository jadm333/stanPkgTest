```R
try(roxygen2::roxygenize(load_code = sourceDir), silent = TRUE)
pkgbuild::compile_dll()
roxygen2::roxygenize()
```
In RStudio press in Build -> Install and Restart
