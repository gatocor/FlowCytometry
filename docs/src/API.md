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
Gating.isInsideGate
Gating.manualGating!
Gating.automaticQC!
```

## Compensation (Spillover/Unmixing)

```@docs
Compensation.computeCompensationMatrix!
Compensation.compensate!
Compensation.assignCompensation!
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