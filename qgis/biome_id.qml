<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34" styleCategories="Symbology">
  <pipe>
    <provider>
      <resampling enabled="false" zoomedInResamplingMethod="nearestNeighbour" zoomedOutResamplingMethod="nearestNeighbour" maxOversampling="2"/>
    </provider>
    <rasterrenderer type="paletted" opacity="1" band="1" alphaBand="-1" nodataColor="">
      <rasterTransparency/>
      <minMaxOrigin>
        <limits>None</limits>
        <extent>WholeRaster</extent>
        <statAccuracy>Estimated</statAccuracy>
        <cumulativeCutLower>0.02</cumulativeCutLower>
        <cumulativeCutUpper>0.98</cumulativeCutUpper>
        <stdDevFactor>2</stdDevFactor>
      </minMaxOrigin>
      <colorPalette>
        <paletteEntry value="0" color="#bcdcee" label="Ocean" alpha="255"/>
        <paletteEntry value="1" color="#f5f7f8" label="Permanent Snow &amp; Ice" alpha="255"/>
        <paletteEntry value="2" color="#8f8579" label="Alpine Fellfield" alpha="255"/>
        <paletteEntry value="3" color="#a89bb0" label="Alpine Tundra" alpha="255"/>
        <paletteEntry value="4" color="#1f6f54" label="Subalpine Wet Forest" alpha="255"/>
        <paletteEntry value="5" color="#4f8a5c" label="Subalpine Woodland" alpha="255"/>
        <paletteEntry value="6" color="#b8622e" label="Subalpine Dry Scrub" alpha="255"/>
        <paletteEntry value="7" color="#2e7d46" label="Temperate Forest" alpha="255"/>
        <paletteEntry value="8" color="#7a9c4a" label="Woodland / Shrubland" alpha="255"/>
        <paletteEntry value="9" color="#e0a83f" label="Lowland Steppe / Grassland" alpha="255"/>
      </colorPalette>
      <colorramp type="randomcolors" name="[source]"/>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <huesaturation colorizeGreen="128" colorizeOn="0" colorizeRed="255" colorizeBlue="128" grayscaleMode="0" saturation="0" colorizeStrength="100" invertColors="0"/>
    <rasterresampler maxOversampling="2"/>
    <resamplingStage>resamplingFilter</resamplingStage>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
