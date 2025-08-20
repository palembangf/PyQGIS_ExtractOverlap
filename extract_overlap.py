"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (QgsProcessing,
                       QgsFeatureSink,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFeatureSink,
                       QgsFeature,
                       QgsFields,
                       QgsField,
                       QgsWkbTypes)
from qgis import processing


class PolygonOverlapAlgorithm(QgsProcessingAlgorithm):
    """
    This algorithm takes a polygon vector layer and creates a new layer
    containing the overlapping areas between polygons.
    
    For each pair of overlapping polygons, it creates a new feature
    representing their intersection geometry, along with the IDs of
    the original overlapping features.
    """

    # Constants used to refer to parameters and outputs
    INPUT = 'INPUT'
    OUTPUT = 'OUTPUT'

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return PolygonOverlapAlgorithm()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'polygonoverlapdetection'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Extract Overlap Features')

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr('Vector analysis')

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'vectoranalysis'

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it.
        """
        return self.tr("Detects overlapping areas between polygons in a vector layer. "
                      "Creates a new layer containing intersection geometries with "
                      "attributes storing the IDs of the overlapping features.")

    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        # We add the input vector features source. It should be polygon geometry.
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                self.tr('Input polygon layer'),
                [QgsProcessing.TypeVectorPolygon]
            )
        )

        # We add a feature sink in which to store our processed features (this
        # usually takes the form of a newly created vector layer when the
        # algorithm is run in QGIS).
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                self.tr('Overlap areas')
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """

        # Retrieve the feature source and sink. The 'dest_id' variable is used
        # to uniquely identify the feature sink, and must be included in the
        # dictionary returned by the processAlgorithm function.
        source = self.parameterAsSource(
            parameters,
            self.INPUT,
            context
        )

        # If source was not found, throw an exception to indicate that the algorithm
        # encountered a fatal error. The exception text can be any string, but in this
        # case we use the pre-built invalidSourceError method to return a standard
        # helper text for when a source cannot be evaluated
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT))

        # Create output fields
        fields = QgsFields()
        fields.append(QgsField('id1', QVariant.Int))
        fields.append(QgsField('id2', QVariant.Int))

        (sink, dest_id) = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            fields,
            QgsWkbTypes.Polygon,
            source.sourceCrs()
        )

        # Send some information to the user
        feedback.pushInfo(f'CRS is {source.sourceCrs().authid()}')

        # If sink was not created, throw an exception to indicate that the algorithm
        # encountered a fatal error. The exception text can be any string, but in this
        # case we use the pre-built invalidSinkError method to return a standard
        # helper text for when a sink cannot be evaluated
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT))

        # Get all features from source
        features = list(source.getFeatures())
        total_comparisons = len(features) * (len(features) - 1) / 2
        
        feedback.pushInfo(f'Processing {len(features)} features, {int(total_comparisons)} comparisons')

        # Compute the number of steps to display within the progress bar
        total = 100.0 / total_comparisons if total_comparisons > 0 else 0
        current_step = 0
        overlap_count = 0

        # Compare each feature with every other feature
        for i in range(len(features)):
            for j in range(i + 1, len(features)):
                # Stop the algorithm if cancel button has been clicked
                if feedback.isCanceled():
                    break
                
                f1 = features[i]
                f2 = features[j]
                geom1 = f1.geometry()
                geom2 = f2.geometry()
                
                if geom1.intersects(geom2):
                    intersection = geom1.intersection(geom2)
                    if not intersection.isEmpty() and intersection.area() > 0:
                        # Create new feature for the overlap
                        new_feat = QgsFeature()
                        new_feat.setGeometry(intersection)
                        new_feat.setAttributes([f1.id(), f2.id()])
                        
                        # Add feature to sink
                        sink.addFeature(new_feat, QgsFeatureSink.FastInsert)
                        overlap_count += 1

                # Update the progress bar
                current_step += 1
                feedback.setProgress(int(current_step * total))
            
            # Break outer loop if cancelled
            if feedback.isCanceled():
                break

        feedback.pushInfo(f'Found {overlap_count} overlapping areas')

        # Return the results of the algorithm
        return {self.OUTPUT: dest_id}