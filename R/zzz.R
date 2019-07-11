
.onLoad <- function(libname, pkgname) {
    base_url <- "http://content.cruk.cam.ac.uk/fmlab/holding2019"
    needsdl_data <- c("counts.datlinger2017_stim.rda", "counts.datlinger2017_unstim.rda", "counts.dixit2016_K562_highmoi.rda", "counts.dixit2016_K562_lowmoi.rda", "counts.jackson2019_CSTARVE.rda", "counts.jackson2019_DIAUXY.rda", "counts.jackson2019_MMD.rda", "counts.jackson2019_MMETOH.rda", "counts.jackson2019_NLIMGLN.rda", "counts.jackson2019_NLIMNH4.rda", "counts.jackson2019_NLIMPRO.rda", "counts.jackson2019_NLIMUREA.rda", "counts.jackson2019_RAPA.rda", "counts.jackson2019_YPD.rda", "counts.jackson2019_YPETOH.rda", "counts.jaitin2016.rda")
    needsdl_dr_crispseq <- c("counts.jaitin2016.mtx.gz")
    needsdl_dr_cropseq <- c("counts.datlinger2017_stim.mtx.gz", "counts.datlinger2017_unstim.mtx.gz")
    needsdl_dr_jackson <- c("counts.jackson2019_CSTARVE.mtx.gz", "counts.jackson2019_DIAUXY.mtx.gz", "counts.jackson2019_MMD.mtx.gz", "counts.jackson2019_MMETOH.mtx.gz", "counts.jackson2019_NLIMGLN.mtx.gz", "counts.jackson2019_NLIMNH4.mtx.gz", "counts.jackson2019_NLIMPRO.mtx.gz", "counts.jackson2019_NLIMUREA.mtx.gz", "counts.jackson2019_RAPA.mtx.gz", "counts.jackson2019_YPD.mtx.gz", "counts.jackson2019_YPETOH.mtx.gz")
    needsdl_dr_perturbseq <- c("counts.dixit2016_K562_highmoi.mtx.gz", "counts.dixit2016_K562_lowmoi.mtx.gz")

    needsdl_dr <- c(needsdl_dr_crispseq, needsdl_dr_cropseq, needsdl_dr_jackson, needsdl_dr_perturbseq)
    needsdl_subdir <- c(rep("crisp-seq", length(needsdl_dr_crispseq)),
                        rep("crop-seq", length(needsdl_dr_cropseq)),
                        rep("jackson2019", length(needsdl_dr_jackson)),
                        rep("perturb-seq", length(needsdl_dr_perturbseq)))

    liblocation <- file.path(libname, pkgname)
    for (file in needsdl_data) {
        path <- file.path(liblocation, "data")
        R.utils::mkdirs(path)
        localfile <- file.path(path, file)
        if (! file.exists(localfile)) {
            utils::download.file(file.path(base_url, "data", file), localfile)
        }
    }
    for (idx in 1:length(needsdl_dr)) {
        file <- needsdl_dr[[idx]]
        path <- file.path(liblocation, "data-raw", needsdl_subdir[[idx]])
        R.utils::mkdirs(path)
        localfile <- file.path(path, file)
        if (! file.exists(localfile)) {
            utils::download.file(file.path(base_url, "data-raw", needsdl_subdir[[idx]], file), localfile)
        }

    }
}
