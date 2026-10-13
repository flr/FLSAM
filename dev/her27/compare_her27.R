#-------------------------------------------------------------------------------
# Compare FLSAM 3.0.0.9005 (old) and 3.1.0 (new) on the 2026 North Sea herring
# assessment (C:/Repositories/2026_her.27.3a47d_assessment, model_sf.R + model_mf.R).
#
# 1. run_version.R runs the assessment + retro with each version in its own R
#    process and makes the diagnostic plots of report_plots.R
# 2. run_extras.R applies the new functionality (new version only)
# 3. this script compares the results and writes output/index.html with the
#    old and new plots side by side
#-------------------------------------------------------------------------------
libRoot   <- "C:/Users/NielsHintzen/Repositories/FLSAM_libs"
assessDir <- "C:/Repositories/2026_her.27.3a47d_assessment"
devDir    <- "C:/Users/NielsHintzen/Repositories/FLSAM/dev/her27"
outRoot   <- file.path(devDir,"output")
nRetro    <- 10
rerun     <- FALSE            # TRUE: rerun both assessments (takes a while)
tol       <- 1e-8

rscript <- file.path(R.home("bin"),"Rscript")
if(rerun | !file.exists(file.path(outRoot,"new","results.rds"))){
  for(v in c("old","new"))
    system2(rscript,c(file.path(devDir,"run_version.R"),file.path(libRoot,ifelse(v=="old","orig","new")),
                      v,outRoot,assessDir,nRetro))
}
if(rerun | !file.exists(file.path(outRoot,"new_extra","extras_log.txt")))
  system2(rscript,c(file.path(devDir,"run_extras.R"),file.path(libRoot,"new"),outRoot,assessDir))

old <- readRDS(file.path(outRoot,"old","results.rds"))
new <- readRDS(file.path(outRoot,"new","results.rds"))

#- Comparison table ----------------------------------------------------------------
numOf  <- function(x){
  if(is.data.frame(x)) x <- x[,sapply(x,is.numeric),drop=FALSE]
  if(is.list(x) && !is.data.frame(x)) return(unlist(lapply(x,numOf)))
  as.numeric(unlist(x))
}
compare <- function(a,b,case,item,note=""){
  if(identical(a,b)) return(data.frame(case=case,item=item,status="identical",max.rel.diff=0,note=note))
  na <- numOf(a); nb <- numOf(b)
  if(length(na)!=length(nb)) return(data.frame(case=case,item=item,status="DIFFERENT (structure)",max.rel.diff=NA,note=note))
  md <- suppressWarnings(max(abs(na-nb)/pmax(abs(nb),1e-12),na.rm=TRUE))
  data.frame(case=case,item=item,status=ifelse(isTRUE(all.equal(na,nb,tolerance=tol)),"equal within tolerance","DIFFERENT"),
             max.rel.diff=signif(md,3),note=note)
}
items <- c("nlogl","nopar","params","stock.n","harvest","residuals","ssb","fbar","rec","catchabilities","obs.var")
tab <- rbind(
  do.call(rbind,lapply(items,function(i) compare(old$sf[[i]],new$sf[[i]],"single fleet",i))),
  do.call(rbind,lapply(items,function(i) compare(old$mf[[i]],new$mf[[i]],"multi fleet",i))),
  do.call(rbind,lapply(names(old$sf.retro),function(y) do.call(rbind,lapply(c("nlogl","ssb","fbar","rec"),function(i)
    compare(old$sf.retro[[y]][[i]],new$sf.retro[[y]][[i]],paste("sf retro",y),i))))),
  do.call(rbind,lapply(names(old$mf.retro),function(y) do.call(rbind,lapply(c("nlogl","ssb","fbar","rec"),function(i)
    compare(old$mf.retro[[y]][[i]],new$mf.retro[[y]][[i]],paste("mf retro",y),i))))),
  do.call(rbind,lapply(c("sf","mf"),function(m) do.call(rbind,lapply(c("ssb","fbar","rec"),function(tp)
    compare(old$mohn[[m]][[tp]]$rho,-new$mohn[[m]][[tp]]$rho,paste(m,"mohns.rho"),paste(tp,"(old vs -new)"),
            "sign reversed in 3.1 (Mohn 1999)")))))
)
write.csv(tab,file.path(outRoot,"comparison_table.csv"),row.names=FALSE)
cat("\nFLSAM",old$version,"vs",new$version,"\n"); print(table(tab$status))
print(subset(tab,status!="identical"),row.names=FALSE)
cat("\nRun times (s):\n"); print(rbind(old=old$timing,new=new$timing))

#- HTML page: old and new plots side by side (images embedded, one portable file) ----
img <- function(f) paste0("data:image/png;base64,",jsonlite::base64_enc(readBin(f,"raw",file.info(f)$size)))
pngs   <- sort(unique(c(sub("\\.png$","",list.files(file.path(outRoot,"old"),"png$")),
                        sub("\\.png$","",list.files(file.path(outRoot,"new"),"png$")))))
extras <- sub("\\.png$","",list.files(file.path(outRoot,"new_extra"),"png$"))
csvTab <- function(df) paste0("<table><tr>",paste0("<th>",names(df),"</th>",collapse=""),"</tr>",
                              paste0(apply(df,1,function(r) paste0("<tr>",paste0("<td>",r,"</td>",collapse=""),"</tr>")),collapse=""),"</table>")
nonId  <- subset(tab,status!="identical")
html <- c("<!doctype html><html><head><meta charset='utf-8'><title>FLSAM old vs new: her.27.3a47d</title>",
"<meta name='viewport' content='width=device-width,initial-scale=1'>",
"<style>body{font-family:Segoe UI,Arial,sans-serif;margin:16px;color:#1d2b33;background:#fff}",
"h1,h2{color:#004868}h2{border-bottom:3px solid #FF800A;padding-bottom:4px;margin-top:36px}",
".pair{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-bottom:28px}.pair img,.one img{width:100%;border:1px solid #cfd8dc}",
".lbl{font-weight:600;color:#004868}.new .lbl{color:#FF800A}h3{color:#00BADC;margin-bottom:6px}",
"table{border-collapse:collapse;font-size:13px}td,th{border:1px solid #cfd8dc;padding:3px 8px}th{background:#004868;color:#fff}",
".ok{color:#31C381;font-weight:700}summary{cursor:pointer;color:#004868;font-weight:600;margin:8px 0}@media(max-width:800px){.pair{grid-template-columns:1fr}}</style></head><body>",
sprintf("<h1>FLSAM %s vs %s &mdash; 2026 North Sea herring assessment</h1>",old$version,new$version),
sprintf("<p>Assessment: <code>%s</code> (model_sf.R and model_mf.R, data_sms2023, addM = 0.02, retro %d years).</p>",assessDir,nRetro),
"<h2>Numerical comparison</h2>",
sprintf("<p><span class='ok'>%d</span> items identical, %d equal within tolerance (%g), %d different.</p>",
        sum(tab$status=="identical"),sum(tab$status=="equal within tolerance"),tol,sum(grepl("DIFFERENT",tab$status))),
if(nrow(nonId)>0) csvTab(nonId) else "<p>All items are bit-identical.</p>",
"<p>Full table: <code>comparison_table.csv</code>.</p>",
"<h2>Diagnostics of report_plots.R: old (left) and new (right)</h2>",
unlist(lapply(pngs,function(p) c(sprintf("<h3>%s</h3><div class='pair'>",p),
  sprintf("<div><div class='lbl'>FLSAM %s</div><img src='%s' alt='%s old'></div>",old$version,img(file.path(outRoot,'old',paste0(p,'.png'))),p),
  sprintf("<div class='new'><div class='lbl'>FLSAM %s</div><img src='%s' alt='%s new'></div></div>",new$version,img(file.path(outRoot,'new',paste0(p,'.png'))),p)))),
"<h2>Residual diagnostics per fleet and age (residual.diagnostics): old (left) and new (right)</h2>",
unlist(lapply(c(sf="single fleet",mf="multi fleet"),function(m){
  pg <- sort(list.files(file.path(outRoot,"old","resdiag"),paste0("^",names(which(c(sf="single fleet",mf="multi fleet")==m)),"_.*png$")))
  c(sprintf("<details><summary>%s: %d pages (click to open)</summary>",m,length(pg)),
    unlist(lapply(pg,function(f) c("<div class='pair'>",
      sprintf("<div><div class='lbl'>FLSAM %s</div><img src='%s'></div>",old$version,img(file.path(outRoot,"old","resdiag",f))),
      sprintf("<div class='new'><div class='lbl'>FLSAM %s</div><img src='%s'></div></div>",new$version,img(file.path(outRoot,"new","resdiag",f)))))),
    "</details>")})),
"<h2>New functionality (FLSAM 3.1 only)</h2>",
"<p>Tables (csv) and the log with forecast, reference points, jitter test and the test of the old saved objects are in <code>new_extra/</code> (<code>extras_log.txt</code>).</p>",
unlist(lapply(extras,function(p) sprintf("<h3>%s</h3><div class='one'><img src='%s' alt='%s'></div>",p,img(file.path(outRoot,'new_extra',paste0(p,'.png'))),p))),
"</body></html>")
writeLines(html,file.path(outRoot,"index.html"))
cat("\nOverview:",file.path(outRoot,"index.html"),"\n")
