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
```

## IO

```@docs
loadFCExperiment
```

## Gating

```@docs
isInsideGate
```

## Dimensionaity reduction

```@docs
DimensionalityReduction.pca!
DimensionalityReduction.umap!
```
## Clustering

```@docs
Clustering.KMeans!
Clustering.KMeansTuning
Clustering.Hierarchical!
```