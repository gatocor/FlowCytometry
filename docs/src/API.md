# API

## Basic Structures

```@docs
FlowCytometryGate
FlowCytometryControl
FlowCytometryExperiment
```

## Operations with structures

```@docs
removeCells
removeCells!
removeChannels
removeChannels!
renameControl!
checkControlNames
```

## IO

```@docs
loadFCExperiment
loadFCControls
```

## Gating

```@docs
isInsideGate
```

## Spilover/Unmixing

```@docs
Compensation.computeSpilloverMatrix!
```

## Dimensionaity reduction

```@docs
DimensionalityReduction.pca!
DimensionalityReduction.umap!
```
## Clustering

```@docs
Clustering.kmeans!
Clustering.kmeansTuning
Clustering.agglomerative!
```