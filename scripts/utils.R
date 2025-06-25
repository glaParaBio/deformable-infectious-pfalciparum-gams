library(data.table)
library(scales)
library(ggplot2)

parse_sample_sheet <- function(RNAseq_master_samples){
    # Clean and reformat RNAseq_master_samples.csv for use in many rules 

    ss <- copy(RNAseq_master_samples)
    # It appears that the value of "Day" is not the same across all samples. We add
    # an index variable to indicate which samples belong to the same subset
    # ss[, day_offset := ifelse(!is.na(deformability), 'deform_exfl', 'only_infect')]


    # From Deepali 13/06/2025
    #  There shouldn't be an exflagellation value listed for any downstream
    #  samples. I see that in this table the up and down samples at each time
    #  point have the same values so somehow someone must have copied them
    #  over. We should go through any other samples where downstream sample
    #  exflagellation values were included and remove them from analysis.
    # ss[, exfl_XA_per_ml := ifelse(deformability == 'down' & !is.na(deformability), NA, exfl_XA_per_ml)]

    # ss <- ss[strain == 'NF54']
    # ss <- ss[!is.na(exfl_XA_per_ml) | !is.na(inf_mosq_percent) | !is.na(oocysts_per_mosq) |!is.na(deformability)]

    ss[, deformability := ifelse(is.na(deformability), NA,
        ifelse(deformability == 'up', 0,
            ifelse(deformability == 'down', 1, -1)))]
    stopifnot(nrow(ss[deformability == -1]) == 0)

    ss[, log2_exfl_XA_per_ml := log2(exfl_XA_per_ml+1)]

    lss <- melt(data= ss[, list(Day, sample_ID, plate, day_offset,
        `log2(Exflagellation/ml)`= log2_exfl_XA_per_ml, 
        `% Infected mosquitoes`= as.numeric(inf_mosq_percent), 
        `Oocysts per mosquito`= as.numeric(oocysts_per_mosq),
        Deformability = deformability)], 
        id.vars= c('sample_ID', 'Day', 'plate', 'day_offset'), variable.name= 'trait', value.name= 'value')
    suppressWarnings(lss[, Day := as.numeric(as.character(Day))])
    lss <- lss[!is.na(value) & !is.na(Day)]
 
    lss[, dicho_value := NA]
    lss[, dicho_value := ifelse(trait == 'Deformability', ifelse(value == 1, 'Deformable', ifelse(value == 0, 'Mixed', NA)), dicho_value)]
    lss[, dicho_value := ifelse(trait == 'log2(Exflagellation/ml)' | trait == 'Oocysts per mosquito', ifelse(value == 0, 'Not detected', 'Detected'), dicho_value)]
    lss[, dicho_value := ifelse(trait == '% Infected mosquitoes', ifelse(value < 50, 'Low infection', 'High infection'), dicho_value)]
    stopifnot(nrow(lss[is.na(dicho_value)]) == 0)
    return(lss)
}

# CREDIT: https://stackoverflow.com/questions/14613355/how-to-get-something-like-matplotlibs-symlog-scale-in-ggplot-or-lattice
symlog_trans <- function(base = 10, thr = 1, scale = 1){
  trans <- function(x)
    ifelse(abs(x) < thr, x, sign(x) * 
             (thr + scale * suppressWarnings(log(sign(x) * x / thr, base))))

  inv <- function(x)
    ifelse(abs(x) < thr, x, sign(x) * 
             base^((sign(x) * x - thr) / scale) * thr)

  breaks <- function(x){
    sgn <- sign(x[which.max(abs(x))])
    if(all(abs(x) < thr))
      pretty_breaks()(x)
    else if(prod(x) >= 0){
      if(min(abs(x)) < thr)
        sgn * unique(c(pretty_breaks()(c(min(abs(x)), thr)),
                       log_breaks(base)(c(max(abs(x)), thr))))
      else
        sgn * log_breaks(base)(sgn * x)
    } else {
      if(min(abs(x)) < thr)
        unique(c(sgn * log_breaks()(c(max(abs(x)), thr)),
                 pretty_breaks()(c(sgn * thr, x[which.min(abs(x))]))))
      else
        unique(c(-log_breaks(base)(c(thr, -x[1])),
                 pretty_breaks()(c(-thr, thr)),
                 log_breaks(base)(c(thr, x[2]))))
    }
  }
  trans_new(paste("symlog", thr, base, scale, sep = "-"), trans, inv, breaks)
}

# CREDIT: https://rdrr.io/github/willgearty/deeptime/f/README.md
#' Transformed and flipped Cartesian coordinate system
#'
#' `coord_trans_flip` behaves similarly to [ggplot2::coord_trans()] in that it
#' occurs after statistical transformation and will affect the visual appearance
#' of geoms. The main difference is that it also flips the x and y coordinates
#' like [ggplot2::coord_flip()].
#'
#' @importFrom ggplot2 ggproto
#' @inheritParams ggplot2::coord_trans
#' @export
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(disp, wt)) +
#'   geom_point() +
#'   coord_trans_flip(x = "log10", y = "log10")
coord_trans_flip <- function(x = "identity", y = "identity",
                             xlim = NULL, ylim = NULL,
                             clip = "on", expand = TRUE) {
  # resolve transformers
  if (is.character(x)) x <- as.trans(x)
  if (is.character(y)) y <- as.trans(y)

  ggproto(NULL, CoordTransFlip,
    trans = list(x = x, y = y),
    limits = list(x = xlim, y = ylim),
    expand = expand,
    clip = clip
  )
}

# copied from ggplot2
flip_axis_labels <- function(x) {
  old_names <- names(x)

  new_names <- old_names
  new_names <- gsub("^x", "z", new_names)
  new_names <- gsub("^y", "x", new_names)
  new_names <- gsub("^z", "y", new_names)

  setNames(x, new_names)
}

#' @rdname coord_trans_flip
#' @format NULL
#' @usage NULL
#' @importFrom ggplot2 ggproto CoordTrans CoordFlip ggproto_parent
#' @export
CoordTransFlip <- ggproto("CoordTransFlip", CoordTrans,
  transform = function(self, data, panel_params) {
    # Need the panel params to be unflipped to correctly transform the data
    panel_params <- flip_axis_labels(panel_params)
    data <- ggproto_parent(CoordTrans, self)$transform(data, panel_params)
    flip_axis_labels(data)
  },
  backtransform_range = function(self, panel_params) {
    un_flipped_range <-
      ggproto_parent(CoordTrans, self)$backtransform_range(panel_params)
    list(x = un_flipped_range$y, y = un_flipped_range$x)
  },
  range = function(self, panel_params) {
    # summarise_layout() expects the original x and y ranges here,
    # not the ones we would get after flipping the axes
    un_flipped_range <- ggproto_parent(CoordTrans, self)$range(panel_params)
    list(x = un_flipped_range$y, y = un_flipped_range$x)
  },
  setup_panel_params = function(self, scale_x, scale_y, params = list()) {
    parent <- ggproto_parent(CoordTrans, self)
    panel_params <- parent$setup_panel_params(scale_x, scale_y, params)
    flip_axis_labels(panel_params)
  },
  labels = function(labels, panel_params) {
    CoordTrans$labels(flip_axis_labels(labels), panel_params)
  },
  setup_layout = function(layout, params) {
    CoordFlip$setup_layout(layout, params)
  },
  modify_scales = function(scales_x, scales_y) {
    CoordFlip$modify_scales(scales_x, scales_y)
  }
)

# CREDIT: https://stats.stackexchange.com/questions/18887/how-to-calculate-a-confidence-interval-for-spearmans-rank-correlation
spearman_CI <- function(x, y, alpha = 0.05, use='pairwise.complete.obs'){
  rs <- cor(x, y, method = "spearman", use)
  n <- sum(complete.cases(x, y))
  sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2)))
}
