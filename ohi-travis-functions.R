calculate_scores <- function(){
  
  # load required libraries
  wd = getwd()
  suppressWarnings(require(ohicore))

  # ensure on draft repo
  library(git2r)
  repo = repository('.')
  checkout(repo, 'draft')
  
  # iterate through all scenarios (by finding layers.csv)
  dirs_scenario = normalizePath(dirname(list.files('.', 'layers.csv', recursive=T, full.names=T)))
  for (dir_scenario in dirs_scenario){ # dir_scenario=dirs_scenario[1]
    
    # set working directory to scenario
    setwd(dir_scenario)
    cat('\n\nCALCULATE SCORES for SCENARIO', basename(dir_scenario), '\n')
    
    # load scenario configuration
    conf <<- Conf('conf')
    
    # run checks on scenario layers
    CheckLayers('layers.csv', 'layers', flds_id=conf$config$layers_id_fields)
    
    # load scenario layers
    layers <<- Layers('layers.csv', 'layers')
    
    # calculate scenario scores
    scores = CalculateAll(conf, layers, debug=F)
    write.csv(scores, 'scores.csv', na='', row.names=F)
    
    # document versions of packages and specifics of ohicore
    cat(
      capture.output(sessionInfo()), '\n\n',
      readLines(file.path(system.file(package='ohicore'), 'DESCRIPTION')), 
      file='session.txt', sep='\n')   
  }  
  
  setwd(wd)
}

create_results <- function(){
  
  library(methods)
  library(ohicore)
  library(tidyr)
  library(dplyr)
  
  # load required libraries
  suppressWarnings(require(ohicore))
  
  # get env
  for (v in c('study_area','default_branch_scenario')){
    if (!exists(v)) assign(v, Sys.getenv(v))
  }
  
  # ensure draft repo
  system('git checkout draft --force')
  
  # iterate through all scenarios (by finding scores.csv)
  wd = getwd() # presumably in top level folder of repo containing scenario folders
  dirs_scenario = normalizePath(dirname(list.files('.', 'scores.csv', recursive=T, full.names=T)))
  for (dir_scenario in dirs_scenario){ # dir_scenario = '~/github/clip-n-ship/alb/alb2014' # dir_scenario = dirs_scenario[1]
    
    # load scenario configuration, layers and scores
    setwd(dir_scenario)
    conf = Conf('conf')
    layers      = Layers('layers.csv', 'layers')
    scores      = read.csv('scores.csv')
    scores_png  = 'reports/figures/scores_400x250.png'
    regions_csv = 'reports/tables/region_titles.csv'
    
    # get goals for flowers, all and specific to weights
    goals.all = arrange(conf$goals, order_color)[['goal']]
    
    # get colors for aster, based on 10 colors, but extended to all goals. subselect for goals.wts
    cols.goals.all = colorRampPalette(RColorBrewer::brewer.pal(10, 'Spectral'), space='Lab')(length(goals.all))
    names(cols.goals.all) = goals.all
    
    # get subgoals and goals, not supragoals, for doing flower plot
    goals_supra = na.omit(unique(conf$goals$parent))
    wts = with(subset(conf$goals, !goal %in% goals_supra, c(goal, weight)), setNames(weight, goal))
    goal_labels = gsub('\\n', '\n', with(conf$goals, setNames(name_flower, goal))[names(wts)], fixed=T)
    
    # region names, ordered by GLOBAL and alphabetical
    rgn_names = rbind(
      data.frame(
        region_id=0,
        rgn_name='GLOBAL',
        rgn_title=study_area),
      SelectLayersData(layers, layers=conf$config$layer_region_labels, narrow=T) %>%
        select(
          region_id=id_num,
          rgn_name=val_chr)  %>%
        mutate(
          rgn_title=rgn_name) %>%
        arrange(rgn_name))
    dir.create(dirname(regions_csv), showWarnings=F, recursive=T)
    write.csv(rgn_names, regions_csv, row.names=F, na='')
    
    # use factors to sort by goal and dimension in scores
    conf$goals = arrange(conf$goals, order_hierarchy)
    scores$goal_label = factor(
      scores$goal,
      levels = c('Index', conf$goals$goal),
      labels = c('Index', ifelse(!is.na(conf$goals$parent),
                                 sprintf('. %s', conf$goals$name),
                                 conf$goals$name)),
      ordered=T)
    scores$dimension_label = factor(
      scores$dimension,
      levels = names(conf$config$dimension_descriptions),
      ordered=T)
    
    # loop through regions
    for (rgn_id in unique(scores$region_id)){ # rgn_id=0
      
      # rgn vars
      rgn_name   = subset(rgn_names, region_id==rgn_id, rgn_name, drop=T)
      flower_png = sprintf('reports/figures/flower_%s.png', gsub(' ','_', rgn_name))
      scores_csv = sprintf('reports/tables/scores_%s.csv', gsub(' ','_', rgn_name))
      
      # create directories, if needed
      dir.create(dirname(flower_png), showWarnings=F)
      dir.create(dirname(scores_csv), showWarnings=F)
      
      # region scores
      g_x = with(subset(scores, dimension=='score' & region_id==rgn_id ),
                 setNames(score, goal))[names(wts)]
      x   = subset(scores, dimension=='score' & region_id==rgn_id & goal == 'Index', score, drop=T)
      
      # flower plot ----
      res=150
      png(flower_png, width=res*7, height=res*7, res=res)
      PlotFlower(
        #main = rgn_name,
        lengths=ifelse(
          is.na(g_x),
          100,
          g_x),
        widths=wts,
        fill.col=ifelse(
          is.na(g_x),
          'grey80',
          cols.goals.all[names(wts)]),
        labels  =ifelse(
          is.na(g_x),
          paste(goal_labels, '-', sep='\n'),
          paste(goal_labels, round(g_x), sep='\n')),
        center=round(x),
        max.length = 100, disk=0.4, label.cex=0.9, label.offset=0.155, cex=2.2, cex.main=2.5)
      dev.off()
      
      if (rgn_id==0){
        res = 72
        png(scores_png, width=400, height=250, res=res)
        par(omi=c(0,0.85,0,0.85))
        PlotFlower(
          lengths=ifelse(
            is.na(g_x),
            100,
            g_x),
          widths=wts,
          fill.col=ifelse(
            is.na(g_x),
            'grey80',
            cols.goals.all[names(wts)]),
          labels = '',
          center=round(x),
          max.length = 100, disk=0.4, label.cex=0.9, label.offset=0.155, cex=2.2, cex.main=2.5)
        dev.off()
        #system(sprintf('open %s', scores_png))
      }
      
      #system(sprintf('convert -density 150x150 %s %s', fig_pdf, fig_png)) # imagemagick's convert
      
      # table csv ----
      scores %>%
        filter(region_id == rgn_id) %>%
        select(goal_label, dimension_label, score) %>%
        spread(dimension_label, score) %>%
        dplyr::rename(' '=goal_label) %>%
        write.csv(scores_csv, row.names=F, na='')
    }
  }
  setwd(wd)
}

create_pages <- function(){
  
  library(yaml)
  library(brew)
  library(ohicore)
  library(dplyr)
  library(knitr)
  library(stringr)
  library(markdown)
  library(rmarkdown)
  library(httr)
  library(git2r)
  merge = base::merge
  diff  = base::diff
        
  # assume in draft branch, get default_branch_scenario set by .travis.yml
  wd = getwd()
  system('git pull; git checkout draft')
  default_branch_scenario  = Sys.getenv('default_branch_scenario')
  study_area               = Sys.getenv('study_area')
  if (default_branch_scenario == '' | study_area == ''){        
    # if not set, then running locally so read in yaml
    travis_yaml = yaml.load_file('.travis.yml')    
    for (var in travis_yaml$env$global){ # var = travis_yaml$env$global[[2]]
      if (is.null(names(var))){
        var_parts = str_trim(str_split(var, '=')[[1]])
        assign(var_parts[1], str_replace_all(var_parts[2], '\"',''))
      }
    }    
    git_owner = 'OHI-Science'
    git_repo  = basename(wd)
    git_slug  = sprintf('%s/%s', git_owner, git_repo)
    git_url   = sprintf('https://github.com/%s', git_slug)  
  } else {
    git_slug  = Sys.getenv('TRAVIS_REPO_SLUG')
    git_owner = str_split(git_slug, '/')[[1]][1]
    git_repo  = str_split(git_slug, '/')[[1]][2]
    git_url   = sprintf('https://github.com/%s', git_slug)    
  }  
  
  # get template brew files
  # update vector: sprintf("'%s'", paste(list.files('~/github/ohi-webapps/results'), collapse="','"))
  dir_brew   = '~/tmp/ohi-webapps'
  unlink(dir_brew, recursive=T)
  dir.create(dir_brew, recursive=T, showWarnings=F)
  for (f in c('navbar.brew.html','regions.brew.md','layers.brew.md','goals.brew.md','scores.brew.md')){
    url_in = file.path('https://raw.githubusercontent.com/OHI-Science/ohi-webapps/master/results', f)
    f_out  = file.path(dir_brew, f)
    writeBin(httr::content(GET(url_in)), f_out)
  }
  
  # clone repo with all branches
  dir_repo = sprintf('~/tmp/%s', git_repo)  
  dir.create(dir_repo, recursive=T, showWarnings=F)
  unlink(dir_repo, recursive=T)
  repo = clone(git_url, normalizePath(dir_repo, mustWork=F))
  setwd(dir_repo)
  
  # archive branches
  dir_archive = sprintf('~/tmp/%s_archive', git_repo)
  dir.create(dir_archive, recursive=T, showWarnings=F)  
  unlink(dir_archive, recursive=T)
  git_branches   = setdiff(sapply(git2r::branches(repo, flags='remote'), function(x) str_replace(x@name, 'origin/', '')), c('HEAD','gh-pages','app'))
  branch_commits = list()
  for (branch in git_branches){ # branch = 'published'        
    checkout(repo, branch, force=T)
    branch_commits[[branch]] = commits(repo)
    dir_branch = file.path(dir_archive, branch)    
    files = list.files(dir_repo, recursive=T)
    for (f in files){ # f = shiny_files[1]
      dir.create(dirname(file.path(dir_branch, f)), showWarnings=F, recursive=T)
      file.copy(file.path(dir_repo, f), file.path(dir_branch, f), overwrite = T, copy.mode=T, copy.date=T) # suppressWarnings)
    }    
  }
  
  # switch to gh-pages branch
  checkout(repo, branch='gh-pages', force=T)
  
  # get list of all branch/scenarios and directory to output
  branch_scenarios = dirname(list.files(dir_archive, 'scores.csv', recursive=T))
  dir_bs_pages = setNames(
    ifelse(
      branch_scenarios == default_branch_scenario, 
      '.', 
      branch_scenarios),
    branch_scenarios)
  
  # iterate over branch/scenarios
  for (branch_scenario in branch_scenarios){ # branch_scenario=branch_scenarios[2]
    
    # get vars
    branch      =  dirname(branch_scenario)
    scenario    = basename(branch_scenario)
    dir_bs_data = normalizePath(file.path(dir_archive, branch_scenario))
    rgns         = file.path(dir_bs_data, 'scores.csv') %>% read.csv %>% select(region_id) %>% unique %>% getElement('region_id')
    layers       = ohicore::Layers(file.path(dir_bs_data, 'layers.csv'), file.path(dir_bs_data, 'layers'))
    conf         = ohicore::Conf(file.path(dir_bs_data, 'conf'))
    # region names, ordered by GLOBAL and alphabetical
    rgns = rbind(
      data.frame(
        id    = 0, 
        name  = 'GLOBAL',
        title = study_area,
        stringsAsFactors=F),
      ohicore::SelectLayersData(layers, layers=conf$config$layer_region_labels, narrow=T) %>%
        select(
          id    = id_num, 
          name  = val_chr) %>%
        mutate(
          title = name)  %>% 
        arrange(title))
    
    # copy results: figures and tables
    dir_data_results  =  file.path(dir_archive, branch_scenario, 'reports')
    dir_pages_results =  file.path('results', branch_scenario)
    dir.create(dir_pages_results, showWarnings=F, recursive=T)
    file.copy(list.files(dir_data_results, full.names=T), dir_pages_results, recursive=T)
    
    # brew markdown files
    for (f_brew in list.files(dir_brew, '.*\\.brew\\.md', full.names=T)){ # f_brew = list.files(dir_brew, '.*\\.brew\\.md', full.names=T)[1]
      section = str_replace(basename(f_brew), fixed('.brew.md'), '')
      branch_scenario_navbar = utils::capture.output({ suppressWarnings(brew(file.path(dir_brew, 'navbar.brew.html'))) })
      f_md = file.path(dir_bs_pages[[branch_scenario]], section, 'index.md')
      dir.create(dirname(f_md), showWarnings=F, recursive=T)
      cat(sprintf('%s -> %s\n', f_brew, f_md))
      suppressWarnings(brew(f_brew, f_md))
    }
    
    # copy regions
    file.copy(file.path(dir_data_results, 'tables/region_titles.csv'), sprintf('_data/regions_%s.csv', str_replace(branch_scenario, '/', '_')), overwrite=T)
  }
  
  # copy scores_400x250.png used in navigation menu
  file.copy(file.path(dir_archive, default_branch_scenario, 'reports/figures/scores_400x250.png'), 'images/scores_400x250.png', overwrite=T)
  
  # push gh-pages
  k = branch_commits[['draft']][[1]]
  system(sprintf('git add -A; git commit -a -m "automatically create_pages from draft commit %0.7s"', k@sha))
  system(sprintf('git push https://${GH_TOKEN}@github.com/%s.git HEAD:gh-pages', git_slug))
  
  # update status ----
  
  # get status repo  depth of 1 only
  if (file.exists('~/tmp/subcountry')){
    system('cd ~/tmp/subcountry; git pull')
  } else {
    dir.create(dirname('tmp'), showWarnings=F, recursive=T)
    system('git clone --depth=1 https://github.com/OHI-Science/subcountry ~/tmp/subcountry')
  }
  csv_status = '~/tmp/subcountry/_data/status.csv'
  d = read.csv(csv_status, stringsAsFactors=F)
  
  # get this repo's info
  n_rgns = file.path(dir_archive, default_branch_scenario, 'reports/tables/region_titles.csv') %>% read.csv() %>% nrow() - 1    
  k = branch_commits[['draft']][[1]]
  
  # update status
  i = which(d$repo == git_repo)
  d$status[i]    = sprintf('[![](https://api.travis-ci.org/OHI-Science/%s.svg?branch=draft)](https://travis-ci.org/OHI-Science/%s/branches)', git_repo, git_repo)
  d$last_mod[i]  = sprintf('%0.10s', as(k@author@when, 'character'))
  d$last_sha[i]  = sprintf('%0.7s', k@sha)
  d$last_msg[i]  = k@summary
  d$map_url[i]   = sprintf('http://ohi-science.org/%s/images/regions_30x20.png', git_repo)
  d$n_regions[i] = n_rgns
 
  # update status repo
  write.csv(d, csv_status, row.names=F, na='')
  system(sprintf('cd ~/tmp/subcountry; git commit -a -m "updated status from %s commit %0.7s"', git_repo, k@sha))
  system('cd ~/tmp/subcountry; git push https://${GH_TOKEN}@github.com/OHI-Science/subcountry.git HEAD:gh-pages')  
  
  # return to original directory
  setwd(wd)
}

push_branch <- function(branch='draft', ci_skip=T){
  
  # set message with [ci skip] to skip travis-ci build for this push
  ci_skip_msg = c('TRUE'='\n[ci skip]', 'FALSE'='')[as.character(ci_skip)]
  
  if (all(Sys.getenv('GH_TOKEN') > '', Sys.getenv('TRAVIS_COMMIT') > '', Sys.getenv('TRAVIS_REPO_SLUG') > '')){
    
    # working on travis-ci
    system(sprintf('git add -A; git commit -a -m "automatically calculate_scores from commit ${TRAVIS_COMMIT}%s"', ci_skip_msg))
    system(sprintf('git push https://${GH_TOKEN}@github.com/${TRAVIS_REPO_SLUG}.git HEAD:%s', branch))
    
  } else {
    
    # working locally, gh_token set in create_init.R, repo_name set in create_init_sc.Rs
    system(sprintf('git add -A; git commit -a -m "automatically calculate_scores from commit `git rev-parse HEAD`%s" ', ci_skip_msg))
    system(sprintf('git push https://%s@github.com/%s.git HEAD:%s', gh_token, git_slug, branch))
    
  }
}

# main
args <- commandArgs(trailingOnly=T)
if (length(args)>0){
  message('')
  fxn = args[1]
  if (length(args)==1){
    eval(parse(text=sprintf('%s()', fxn)))
  } else {
    eval(parse(text=sprintf("%s('%s')", fxn, paste( args[2:length(args)], collapse="', '"))))
  }
}