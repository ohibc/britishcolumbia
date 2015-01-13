# load required libraries
suppressWarnings(require(ohicore))

# set working directory to the scenario directory, ie containing conf and layers directories
setwd('/home/bbest/github/clip-n-ship/can/subcountry2014')

# load scenario configuration
conf = Conf('conf')

# run checks on scenario layers
CheckLayers('layers.csv', 'layers', flds_id=conf$config$layers_id_fields)

# load scenario layers
layers = Layers('layers.csv', 'layers')

# calculate scenario scores
scores = CalculateAll(conf, layers, debug=F)
write.csv(scores, 'scores.csv', na='', row.names=F)

## testing this
