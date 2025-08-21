rm(list=ls())
setwd('~/Dropbox (University of Oregon)/')
setwd("parasite_networks")

## Script for plotting all of the important explanatory variables.
library(tidyverse)
library(viridis)
library(networkD3)
library(htmlwidgets)
library(webshot2)
library(magick)

source("src/misc.R")
source("src/ggplotThemes.R")

load("../parasite_networks/data/allNets.RData")
load(file="../parasite_networks/data/sp_mets.RData")



sp.network.metrics$NetworkKey <- paste(sp.network.metrics$Site, sp.network.metrics$Year,
                             sp.network.metrics$SampleRound, sep=".")

length(all.nets)
all.nets <- all.nets[names(all.nets) %in% sp.network.metrics$NetworkKey]
length(all.nets)

# ---- Utilities ----
scale01 <- function(v) {
  v <- as.numeric(v)
  if (!length(v) || all(!is.finite(v))) return(rep(0.5, length(v)))
  r <- range(v, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || diff(r) < 1e-12) return(rep(0.5, length(v)))
  (v - r[1]) / (r[2] - r[1])
}

parse_net_name <- function(x) {
  parts <- strsplit(x, "\\.")[[1]]
  if (length(parts) < 3) stop(sprintf("Name '%s' not 'Site.Year.SampleRound'.", x))
  list(Site = trimws(parts[1]),
       Year = suppressWarnings(as.integer(parts[2])),
       SampleRound = suppressWarnings(as.integer(parts[3])))
}

clean_sp_metrics <- function(sp.network.metrics) {
  sp.network.metrics %>%
    mutate(
      Site = trimws(as.character(Site)),
      Year = suppressWarnings(as.integer(Year)),
      SampleRound = suppressWarnings(as.integer(SampleRound)),
      GenusSpecies = trimws(as.character(GenusSpecies))
    )
}

get_event_slice <- function(spm, site, year, round) {
  spm %>%
    filter(
      Site == trimws(as.character(site)),
      Year == suppressWarnings(as.integer(year)),
      SampleRound == suppressWarnings(as.integer(round))
    )
}

# Get species -> numeric vector from a given column, averaged within the event
species_vector_from_col <- function(event_df, colname, species) {
  if (is.null(colname)) return(setNames(rep(NA_real_, length(species)), species))
  if (!colname %in% names(event_df)) {
    warning(sprintf("Column '%s' not found in sp.network.metrics; returning NA.", colname))
    return(setNames(rep(NA_real_, length(species)), species))
  }
  df <- event_df %>%
    mutate(.val = suppressWarnings(as.numeric(.data[[colname]]))) %>%
    group_by(GenusSpecies) %>%
    summarize(val = mean(.val, na.rm = TRUE), .groups = "drop")
  out <- df$val[match(species, df$GenusSpecies)]
  setNames(as.numeric(out), species)
}

# ---- Prevalence & Pollinator Size ----
get_prevalence <- function(sp.network.metrics, site, year, round, species, prevalence_col) {
  spm <- clean_sp_metrics(sp.network.metrics)
  ev  <- get_event_slice(spm, site, year, round)
  species <- trimws(as.character(species))
  if (!nrow(ev)) return(setNames(rep(NA_real_, length(species)), species))
  species_vector_from_col(ev, prevalence_col, species)
}

get_pollinator_size <- function(sp.network.metrics, site, year, round, species, pollinator_size_col) {
  spm <- clean_sp_metrics(sp.network.metrics)
  ev  <- get_event_slice(spm, site, year, round)
  species <- trimws(as.character(species))
  if (!nrow(ev)) return(setNames(rep(NA_real_, length(species)), species))
  species_vector_from_col(ev, pollinator_size_col, species)
}

# ---- Matrix -> nodes/links; compute plant/poll sizes; attach labels ----
bipartite_to_d3 <- function(mat,
                            prev_vec = NULL,
                            show_labels = TRUE,
                            plant_size_mode = c("strength","degree"),
                            poll_size_vec = NULL,
                            poll_fallback_mode = c("strength","degree")) {
  plant_size_mode    <- match.arg(plant_size_mode)
  poll_fallback_mode <- match.arg(poll_fallback_mode)

  if (inherits(mat, "data.frame")) mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  nR <- nrow(mat); nC <- ncol(mat)

  if (is.null(rownames(mat))) rownames(mat) <- paste0("Plant_", seq_len(nR))
  if (is.null(colnames(mat))) colnames(mat) <- paste0("Poll_",  seq_len(nC))

  plants <- rownames(mat)
  polls  <- colnames(mat)

  # Plant metrics
  plant_strength <- if (nR) rowSums(mat, na.rm = TRUE) else numeric(0)
  plant_degree   <- if (nR) rowSums(mat != 0) else numeric(0)
  plant_size_raw <- switch(plant_size_mode,
                           strength = plant_strength,
                           degree   = plant_degree)
  plant_size <- scale01(plant_size_raw)

  # Pollinator metrics
  poll_strength <- if (nC) colSums(mat, na.rm = TRUE) else numeric(0)
  poll_degree   <- if (nC) colSums(mat != 0) else numeric(0)

  if (!is.null(poll_size_vec) && length(poll_size_vec)) {
    hit <- match(polls, names(poll_size_vec))
    poll_size_raw <- poll_size_vec[hit]
  } else {
    poll_size_raw <- switch(poll_fallback_mode,
                            strength = poll_strength,
                            degree   = poll_degree)
  }
  if (all(is.na(poll_size_raw))) {
    poll_size_raw <- switch(poll_fallback_mode,
                            strength = poll_strength,
                            degree   = poll_degree)
  }
  poll_size <- scale01(poll_size_raw)

  # Nodes
  nodes <- tibble(
    name     = c(plants, polls),
    role     = c(rep("plant", length(plants)), rep("pollinator", length(polls))),
    label    = if (show_labels) c(plants, polls) else rep("", length(plants) + length(polls)),
    prev     = NA_real_,
    pstr     = c(scale01(plant_strength), rep(NA_real_, length(polls))),  # for green shades
    nodesize = c(plant_size, poll_size)
  )

  # Fill pollinator prevalence by NAME
  if (!is.null(prev_vec) && length(prev_vec)) {
    hit <- match(polls, names(prev_vec))
    nodes$prev[(length(plants)+1):(length(plants)+length(polls))] <- prev_vec[hit]
  }

  # Links
  if (nR == 0 || nC == 0) {
    links <- tibble(source = integer(0), target = integer(0), value = numeric(0))
  } else {
    nz <- which(mat != 0, arr.ind = TRUE)
    if (length(nz) == 0) {
      links <- tibble(source = integer(0), target = integer(0), value = numeric(0))
    } else {
      r <- nz[,1]; c <- nz[,2]
      links <- tibble(
        source = as.integer(r - 1L),
        target = as.integer(nR + c - 1L),
        value  = as.numeric(mat[cbind(r,c)])
      )
    }
  }

  list(nodes = nodes, links = links)
}

# ---- Colors: plants -> green; pollinators -> magma(0..1); missing -> gray ----
make_colour_groups_magma <- function(nodes) {
  groups <- rep("NA", nrow(nodes))

  # Plants colored by pstr (row strength) shades of green
  is_plant <- nodes$role == "plant"
  if (any(is_plant)) {
    pstr <- pmin(pmax(nodes$pstr[is_plant], 0), 1)
    groups[is_plant] <- sprintf("PLANT_%03d", floor(pstr * 100))
  }

  # Pollinators with prevalence use magma
  is_poll_prev <- nodes$role == "pollinator" & !is.na(nodes$prev)
  if (any(is_poll_prev)) {
    prev <- pmin(pmax(nodes$prev[is_poll_prev], 0), 1)
    groups[is_poll_prev] <- sprintf("PREV_%03d", floor(prev * 100))
  }

  nodes$group <- groups

  greens <- grDevices::colorRampPalette(c("#e5f5e0", "#74c476", "#006d2c"))(101)  # 0..100
  magma  <- viridisLite::magma(101)

  domain <- c(
    paste0("\"PLANT_", sprintf("%03d", 0:100), "\""),
    "\"NA\"",
    paste0("\"PREV_",  sprintf("%03d", 0:100), "\"")
  )
  range <- c(
    paste0("\"", greens, "\""),
    "\"#bdbdbd\"",
    paste0("\"", magma,  "\"")
  )

  colourScale <- htmlwidgets::JS(
    sprintf("d3.scaleOrdinal().domain([%s]).range([%s])",
            paste(domain, collapse=","), paste(range, collapse=",")) )

  list(nodes = nodes, colourScale = colourScale, magma = magma)
}

# ---- Render ONE widget (labels + prevalence colorbar + size legend) ----
render_one_network <- function(mat, nm, prev_vec,
                               show_labels = TRUE,
                               plant_size_mode = c("strength","degree"),
                               poll_size_vec = NULL,
                               poll_fallback_mode = c("strength","degree"),
                               legend_title = "Parasite prevalence",
                               size_legend_title = "Node size (scaled)") {

  # These MUST match the forceNetwork radiusCalculation constants
  baseRadius  <- 4
  scaleRadius <- 10

  d3dat   <- bipartite_to_d3(mat,
                             prev_vec = prev_vec,
                             show_labels = show_labels,
                             plant_size_mode = plant_size_mode,
                             poll_size_vec = poll_size_vec,
                             poll_fallback_mode = poll_fallback_mode)
  colinfo <- make_colour_groups_magma(d3dat$nodes)
  d3dat$nodes <- colinfo$nodes

  nodeCharge <- htmlwidgets::JS("function(d){return d.role==='plant'?-80:-120;}")

  widg <- networkD3::forceNetwork(
    Links = d3dat$links,
    Nodes = d3dat$nodes,
    Source = "source", Target = "target", Value = "value",
    NodeID = "name",                    # tooltip (labels added below)
    Group  = "group",
    Nodesize = "nodesize",
    radiusCalculation = htmlwidgets::JS(sprintf("%d + %d*d.nodesize", baseRadius, scaleRadius)),
    opacity = 0.9, zoom = TRUE, bounded = TRUE,
    linkDistance = 30, charge = nodeCharge, fontSize = 12,
    colourScale = colinfo$colourScale
  )

                                        # 1) Persistent labels (optional)
  if (show_labels) {
    widg <- htmlwidgets::onRender(
      widg,
      "
      function(el, x) {
        var svg = d3.select(el).select('svg');

        function nodeSel() {
          var n = svg.selectAll('circle.node');
          if (n.empty()) n = svg.selectAll('g.node');
          return n;
        }

        var labelsLayer = svg.select('g.labels');
        if (labelsLayer.empty()) labelsLayer = svg.append('g').attr('class', 'labels');

        function updateLabels() {
          var nodes = nodeSel();
          var data  = nodes.data();
          var t = labelsLayer.selectAll('text.nodelabel').data(data);

          t.enter().append('text')
            .attr('class','nodelabel')
            .style('font-size','10px')
            .style('fill','#111')
            .style('pointer-events','none')
            .merge(t)
            .text(function(d){ return d.label || d.name || ''; })
            .attr('x', function(d){ return (d.x || 0) + 8; })
            .attr('y', function(d){ return (d.y || 0) + 3; });

          t.exit().remove();
          labelsLayer.raise();
        }

        var count = 0, maxCount = 30;
        var iv = setInterval(function(){
          updateLabels();
          if (++count >= maxCount) clearInterval(iv);
        }, 80);

        setTimeout(updateLabels, 600);
      }
      "
    )
  }
# --- Replace your legend onRender(...) with this dynamic, network-anchored version ---
magma_vec  <- colinfo$magma
magma_js   <- paste0("[", paste(sprintf('\"%s\"', magma_vec), collapse=","), "]")
title_js   <- jsonlite::toJSON(legend_title, auto_unbox = TRUE)
size_title <- jsonlite::toJSON(size_legend_title, auto_unbox = TRUE)
base_js    <- jsonlite::toJSON(baseRadius, auto_unbox = TRUE)
scale_js   <- jsonlite::toJSON(scaleRadius, auto_unbox = TRUE)

widg <- htmlwidgets::onRender(
  widg,
  sprintf("
  function(el, x){
    var svg = d3.select(el).select('svg');

    function nodeSel(){
      var n = svg.selectAll('circle.node');
      if (n.empty()) n = svg.selectAll('g.node');
      return n;
    }

    // Draw legends near the node cloud (after layout has mostly settled)
    function placeLegend(){
      var nodes = nodeSel();
      var data  = nodes.data();
      if (!data || !data.length) return;

      // Node bounding box
      var minX = d3.min(data, function(d){ return d.x; }),
          maxX = d3.max(data, function(d){ return d.x; }),
          minY = d3.min(data, function(d){ return d.y; }),
          maxY = d3.max(data, function(d){ return d.y; });

      // SVG size (fallback to DOM box if attrs missing)
      var box = el.getBoundingClientRect();
      var svgW = +svg.attr('width')  || box.width;
      var svgH = +svg.attr('height') || box.height;

      // Legend layout (vertical colorbar + vertical size legend)
      var barW = 12, barH = 180, gap = 18, pad = 12;

      // Size legend height estimate for 3 circles (0, 0.5, 1)
      var baseR = %s, scaleR = %s;
      var svals = [0, 0.5, 1.0];
      var radii = svals.map(function(s){ return baseR + scaleR*s; });
      var sizeLegendH = radii.reduce(function(acc, r){ return acc + (r*2); }, 0) + (svals.length-1)*10 + 22;

      var legendW = Math.max(barW + 40, 110); // width incl. tick labels
      var legendH = barH + gap + sizeLegendH;

      // Try left side of node cloud first; if clipped, use right side
      var left = minX - legendW - pad;
      if (left < pad) left = Math.min(svgW - legendW - pad, maxX + pad);
      var top  = minY;  // top aligned to top of node cloud

      // Clamp inside the SVG
      left = Math.max(pad, Math.min(left, svgW - legendW - pad));
      top  = Math.max(pad, Math.min(top,  svgH - legendH - pad));

      // Clear old legends and draw fresh
      svg.selectAll('g.legend-root, g.legend, defs #magmaGradV2').remove();

      var root = svg.append('g').attr('class','legend-root')
                    .attr('transform','translate('+left+','+top+')');

      // ----- Vertical prevalence colorbar -----
      var defs = svg.select('defs'); if (defs.empty()) defs = svg.append('defs');
      var grad = defs.append('linearGradient')
                     .attr('id','magmaGradV2')
                     .attr('x1','0%%').attr('y1','100%%')
                     .attr('x2','0%%').attr('y2','0%%');
      var pal = %s;  // magma palette
      for (var i=0;i<pal.length;i++){
        grad.append('stop')
            .attr('offset', (i/(pal.length-1)*100)+'%%')
            .attr('stop-color', pal[i]);
      }

      var g1 = root.append('g').attr('class','legend legend-prevalence');
      g1.append('rect')
        .attr('width',  barW)
        .attr('height', barH)
        .style('fill','url(#magmaGradV2)');

      // tick labels: 0 (bottom) and 1 (top)
      g1.append('text')
        .attr('x', barW + 6).attr('y', barH)
        .attr('dominant-baseline','middle')
        .style('font-size','10px').text('0');
      g1.append('text')
        .attr('x', barW + 6).attr('y', 0)
        .attr('dominant-baseline','middle')
        .style('font-size','10px').text('1');

      g1.append('text')
        .attr('x', barW/2).attr('y', -6)
        .attr('text-anchor','middle')
        .style('font-size','10px')
        .text(%s);

      // ----- Vertical size legend (three gray circles) -----
      var g2 = root.append('g').attr('class','legend legend-size')
                   .attr('transform','translate(0,'+(barH + gap)+')');

      g2.append('text')
        .attr('x', barW/2).attr('y', -6)
        .attr('text-anchor','middle')
        .style('font-size','10px')
        .text(%s);

      var yCursor = 8, centerX = barW/2;
      for (var i=0;i<radii.length;i++){
        var r  = radii[i];
        var cy = yCursor + r;

        g2.append('circle')
          .attr('cx', centerX).attr('cy', cy).attr('r', r)
          .style('fill','#bdbdbd').style('stroke','#888').style('stroke-width',1);

        g2.append('text')
          .attr('x', centerX + r + 8).attr('y', cy)
          .attr('dominant-baseline','middle')
          .style('font-size','10px')
          .text(svals[i]);

        yCursor += (r*2) + 10;
      }
    }

    // Place once after layout warms up, and once more later to be safe
    setTimeout(placeLegend, 900);
    setTimeout(placeLegend, 1600);
  }
  ", base_js, scale_js, magma_js, title_js, size_title)
)

  widg
}

# ---- Main: render all networks to two-per-page PDF ----
render_all_networks_pdf <- function(all.nets, sp.network.metrics,
                                    out_pdf              = "all_networks_forceNetworks.pdf",
                                    prevalence_col       = "PropSpCrithidiaPresence",
                                    pollinator_size_col  = NULL,   # e.g., a metric col in sp.network.metrics
                                    plant_size           = c("strength","degree"),
                                    pollinator_fallback  = c("strength","degree"),
                                    show_labels          = TRUE,
                                    vwidth = 900, vheight = 700) {

  plant_size          <- match.arg(plant_size)
  pollinator_fallback <- match.arg(pollinator_fallback)

  out_dir <- file.path(tempdir(), paste0("forceNets_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  errors <- list()

  png_files <- purrr::imap_chr(all.nets, function(mat, nm) {
    tryCatch({
      mat <- as.matrix(mat)
      colnames(mat) <- trimws(colnames(mat))
      rownames(mat) <- trimws(rownames(mat))

      key <- parse_net_name(nm)
      species  <- colnames(mat)

      # prevalence coloring
      prev_vec <- get_prevalence(sp.network.metrics, key$Site, key$Year, key$SampleRound,
                                 species, prevalence_col = prevalence_col)

      # pollinator node sizes: from sp.network.metrics column (if provided), else network fallback
      poll_size_vec <- NULL
      if (!is.null(pollinator_size_col)) {
        poll_size_vec <- get_pollinator_size(sp.network.metrics, key$Site, key$Year, key$SampleRound,
                                             species, pollinator_size_col = pollinator_size_col)
      }

      # Legend titles
      legend_title     <- paste0(prevalence_col, " (0–1)")
      size_legend_title <- if (!is.null(pollinator_size_col)) {
        paste0(pollinator_size_col, " (scaled)")
      } else {
        paste0("Pollinator ", pollinator_fallback, " (scaled)")
      }

      widg <- render_one_network(mat, nm, prev_vec,
                                 show_labels = show_labels,
                                 plant_size_mode = plant_size,
                                 poll_size_vec = poll_size_vec,
                                 poll_fallback_mode = pollinator_fallback,
                                 legend_title = legend_title,
                                 size_legend_title = size_legend_title)

      safe_nm  <- stringr::str_replace_all(nm, "[^A-Za-z0-9_\\-\\.]", "_")
      html_file <- file.path(out_dir, paste0(safe_nm, ".html"))
      png_file  <- file.path(out_dir, paste0(safe_nm, ".png"))

      htmlwidgets::saveWidget(widg, file = html_file, selfcontained = TRUE)
      # Delay to ensure labels & legends render before snapshot
      webshot2::webshot(html_file, file = png_file, vwidth = vwidth, vheight = vheight,
                        cliprect = "viewport", delay = 2.6)

      # Title banner
      img <- magick::image_read(png_file)
      banner <- magick::image_blank(width = magick::image_info(img)$width, height = 80, color = "white")
      banner <- magick::image_annotate(banner, nm, gravity = "center", size = 32, weight = 600)
      titled <- magick::image_append(c(banner, img), stack = TRUE)
      magick::image_write(titled, path = png_file, format = "png")

      png_file
    }, error = function(e) {
      message("Skipping '", nm, "' due to: ", conditionMessage(e))
      errors[[nm]] <<- conditionMessage(e)
      NA_character_
    })
  })

  png_files <- png_files[!is.na(png_files)]
  if (!length(png_files)) stop("All networks failed to render. See errors printed above.")

  # two per page
  imgs <- magick::image_read(png_files)
  pages <- list(); i <- 1
  while (i <= length(imgs)) {
    top <- imgs[i]
    if (i + 1 <= length(imgs)) {
      bottom <- imgs[i + 1]
    } else {
      info <- magick::image_info(top)
      bottom <- magick::image_blank(width = info$width, height = info$height, color = "white")
    }
    pages[[length(pages)+1]] <- magick::image_append(c(top, bottom), stack = TRUE)
    i <- i + 2
  }
  pdf_stack <- do.call(magick::image_join, pages)
  magick::image_write(pdf_stack, path = out_pdf, format = "pdf")

  if (length(errors)) message("Some networks were skipped: ", paste(names(errors), collapse = ", "))
  message("Wrote PDF: ", normalizePath(out_pdf))
  invisible(out_pdf)
}


out_file <- render_all_networks_pdf(
  all.nets, sp.network.metrics,
  out_pdf             = "figures/networks/PropSpCrithidiaPresence_Closeness_forceNetworks.pdf",
  prevalence_col      = "PropSpCrithidiaPresence",   # or "PropSpApicystisSpp"
  pollinator_size_col = "zweighted.closeness",                        # or a metric col from sp.network.metrics
  plant_size          = "strength",                  # "strength" or "degree"
  pollinator_fallback = "degree",                    # "strength" or "degree"
  show_labels         = TRUE
)


out_file <- render_all_networks_pdf(
  all.nets, sp.network.metrics,
  out_pdf             = "figures/networks/PropSpApicystisSpp_Closeness_forceNetworks.pdf",
  prevalence_col      = "PropSpApicystisSpp",   
  pollinator_size_col = "zweighted.closeness",                        # or a metric col from sp.network.metrics
  plant_size          = "strength",                  # "strength" or "degree"
  pollinator_fallback = "degree",                    # "strength" or "degree"
  show_labels         = TRUE
)



out_file <- render_all_networks_pdf(
  all.nets, sp.network.metrics,
  out_pdf             = "figures/networks/PropSpCrithidiaPresence_Degree_forceNetworks.pdf",
  prevalence_col      = "PropSpCrithidiaPresence",   # or "PropSpApicystisSpp"
  pollinator_size_col = "zdegree",                        # or a metric col from sp.network.metrics
  plant_size          = "strength",                  # "strength" or "degree"
  pollinator_fallback = "degree",                    # "strength" or "degree"
  show_labels         = TRUE
)


out_file <- render_all_networks_pdf(
  all.nets, sp.network.metrics,
  out_pdf             = "figures/networks/PropSpApicystisSpp_Degree_forceNetworks.pdf",
  prevalence_col      = "PropSpApicystisSpp",   
  pollinator_size_col = "zdegree",                        # or a metric col from sp.network.metrics
  plant_size          = "strength",                  # "strength" or "degree"
  pollinator_fallback = "degree",                    # "strength" or "degree"
  show_labels         = TRUE
)
