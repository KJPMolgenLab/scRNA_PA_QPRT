library(workflowr)

#wflow_start("./", existing = TRUE)

options(workflowr.sysgit = "")
wflow_build("./analysis/*.Rmd", update=T)
wflow_publish(c("./analysis/*.Rmd", "./code/*", "./docs/*"))
system("git push -u origin main")

