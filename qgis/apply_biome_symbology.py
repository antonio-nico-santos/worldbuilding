"""
Run inside the QGIS Python Console (Plugins -> Python Console) to apply
Tappa 5's biome palette programmatically -- equivalent to loading
qgis/biome_id.qml via Layer Properties -> Symbology -> Style -> Load Style,
but scripted so it stays reproducible and easy to re-apply if the palette
in src/biomes/world_biomes.py (BIOME_NAMES / BIOME_COLORS_HEX) ever changes.
Values/colors/labels below are copied from that module -- keep them in sync
if you edit one, not the other.
"""

from qgis.core import QgsRasterLayer, QgsPalettedRasterRenderer, QgsProject
from qgis.PyQt.QtGui import QColor

# id -> (name, hex) -- see docs/decisions/05_tappa5_biomes.md S5 for the
# belt x moisture derivation of each class.
BIOME_NAMES = [
    "Ocean",
    "Permanent Snow & Ice",
    "Alpine Fellfield",
    "Alpine Tundra",
    "Subalpine Wet Forest",
    "Subalpine Woodland",
    "Subalpine Dry Scrub",
    "Temperate Forest",
    "Woodland / Shrubland",
    "Lowland Steppe / Grassland",
]
BIOME_COLORS_HEX = [
    "#bcdcee",
    "#f5f7f8",
    "#8f8579",
    "#a89bb0",
    "#1f6f54",
    "#4f8a5c",
    "#b8622e",
    "#2e7d46",
    "#7a9c4a",
    "#e0a83f",
]

# Path used if the layer isn't already loaded in this project. Adjust if
# your project keeps the raster elsewhere.
RASTER_PATH = r"C:\projects\worldbuilding\data\processed\biomes\biome_id.bin"
LAYER_NAME = "biome_id"


def get_or_load_layer() -> QgsRasterLayer:
    existing = QgsProject.instance().mapLayersByName(LAYER_NAME)
    if existing:
        return existing[0]
    layer = QgsRasterLayer(RASTER_PATH, LAYER_NAME)
    if not layer.isValid():
        raise RuntimeError(f"Could not load raster at {RASTER_PATH} -- check the path.")
    # Same CRS caveat as every Tappa: the .hdr's map info supplies the affine
    # georeferencing but not a full CRS. Assign "Fictional World LCC" from
    # Layer Properties -> Source -> Assigned CRS if QGIS loads it unknown.
    QgsProject.instance().addMapLayer(layer)
    return layer


def apply_biome_palette(layer: QgsRasterLayer) -> None:
    classes = [
        QgsPalettedRasterRenderer.Class(i, QColor(hexcolor), name)
        for i, (name, hexcolor) in enumerate(zip(BIOME_NAMES, BIOME_COLORS_HEX))
    ]
    renderer = QgsPalettedRasterRenderer(layer.dataProvider(), 1, classes)
    layer.setRenderer(renderer)
    layer.triggerRepaint()


layer = get_or_load_layer()
apply_biome_palette(layer)
print(f"Applied Tappa 5 biome palette ({len(BIOME_NAMES)} classes) to '{layer.name()}'.")
