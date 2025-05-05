// Cargar Sentinel 1 VV
var s1vv = ee.ImageCollection('COPERNICUS/S1_GRD')
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
                    .filter(ee.Filter.eq('instrumentMode', 'IW'))
                    .select('VV')
                    .filterDate('2020-05-10', '2020-05-11')
                    .filterBounds(aoi)
                    .map(function(image) {
                      var edge = image.lt(-30.0);
                      var maskedImage = image.mask().and(edge.not());
                      return image.updateMask(maskedImage);
        });
//Filtro VV por tipo de órbita
var desc = s1vv.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));
var asc = s1vv.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));

print(s1vv,"s1vv")


// Clasifica las imágenes por fecha (más recientes primero) 
var sortedCollection = s1vv.sort('system:time_start', false);

// Selecciona la imagen más reciente que cubre el área de estudio 
var bestImage = sortedCollection.first()
                             .clip(aoi);
print(bestImage,"vv")
// Muestra la imagen seleccionada en el mapa
Map.centerObject(aoi); 
Map.addLayer(bestImage, {min: -30, max: 10}, 'VV');

// Cargar Sentinel 1 VH
var s1vh = ee.ImageCollection('COPERNICUS/S1_GRD')
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
                    .filter(ee.Filter.eq('instrumentMode', 'IW'))
                    .select('VH')
                    .filterDate('2020-05-10', '2020-05-11')
                    .filterBounds(aoi)
                    .map(function(image) {
                      var edge = image.lt(-30.0);
                      var maskedImage = image.mask().and(edge.not());
                      return image.updateMask(maskedImage);
        });
//Filtro VH por tipo de órbita
var desc = s1vh.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));
var asc = s1vh.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));

// Clasifica las imágenes por fecha (más recientes primero) 
var vh_class = s1vh.sort('system:time_start', false);

// Selecciona la imagen más reciente que cubre el área de estudio 
var vh_new = vh_class.first()
                             .clip(aoi);
print(vh_new,"vh")
// Muestra la imagen seleccionada en el mapa
Map.centerObject(aoi); 
Map.addLayer(vh_new, {min: -30, max: 10}, "VH");


//* Function to mask clouds using the Sentinel-2 QA band
// * @param {ee.Image} image Sentinel-2 image
// * @return {ee.Image} cloud masked Sentinel-2 image
// */

function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000);
}

// Cargar Sentinel 2
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filterDate('2020-05-20', '2020-05-21')
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20))
                  .map(maskS2clouds)
                  .select("B2","B3","B4","B8","B11")
                  .filterBounds(aoi);
print(s2,"s2")
var visualization = {
  min: 0.0,
  max: 0.3,
  bands: ['B4', 'B3', 'B2'],
};

//Seleccionar imagen a analizar
var img = s2.mosaic();

//Resample SWIR1 a 10 m
var b11 = img.select("B11")
var resampled_b11 = b11.resample("bilinear")
                       .reproject({
                         crs:b11.projection(),
                         scale:10
                       })
var resolution = resampled_b11.projection().nominalScale();
print(resolution,"b11 resapling")

Map.addLayer(img, visualization, 'RGB');

//Añadir índices
//El Combined Mangrove Recognition Index (CMRI) es un índice que se utiliza en la detección de manglares en imágenes satelitales
var cmri = img.expression('((NIR - R) / (NIR + R)) - ((G - NIR) / (G + NIR))', {
    'NIR': img.select('B8'),
    'R': img.select('B4'),
    'G': img.select('B3')
}).rename('CMRI');
  // MVI Mangrove Vegetation Index
  var mvi = img.expression('(NIR - GREEN) / (SWIR1 - GREEN)', {
    'NIR': img.select('B8'),
    'GREEN': img.select('B3'),
    'SWIR1': resampled_b11
  }).rename('MVI');
   // MDI Mangrove damage index
  var mdi = img.expression('((NIR - SWIR1) / (NIR * SWIR1)) * 10000', {
    'NIR': img.select('B8'),
    'SWIR1': resampled_b11
  }).rename('MDI');
  // NDWI 
  var ndwi = img.expression('((G - NIR) / (G * NIR))', {
    'G': img.select('B3'),
    'NIR': img.select('B8')
  }).rename('NDWI');
  // SAVI 
  var savi = img.expression('((NIR – RED) / (NIR + RED + 0.5)) * (1 + 0.5)', {
    'NIR': img.select('B8'),
    'RED': img.select('B4')
  }).rename('SAVI');
  // EVI
//  var evi = img.expression('2.5 * ((NIR – RED) / ((NIR) + (6 * RED) – (7.5 * BLUE) + 1))', {
//    'NIR': img.select('B8'),
//    'RED': img.select('B4'),
//    "BLUE":img.select('B2')
//  }).rename('EVI');
// Calcular el Índice de Vegetación de Diferencia Normalizada (NDVI)
var ndvi = img.normalizedDifference(['B8','B4']).rename('NDVI');
// Calcular el Índice de Agua Normalizada (NDMI)
//var ndmi = img.normalizedDifference(['B8',resampled_b11]).rename('NDMI');

//Cargar DEM
var dem = ee.Image('USGS/SRTMGL1_003');
var elevation = dem.select('elevation').clip(aoi);


//Cargar imagen con índices
var s2_index = img.addBands([ndvi, cmri, mvi, mdi, ndwi,elevation, bestImage, vh_new]);

print(s2_index,"s2_index")

// Segmentation -----------------------------------------------------------------------------

var sample = s2_index.sample({
  region : aoi,
  scale :30,
  numPixels: 5000
});

//Configurar parametros de algoritmo de K Means

var numClusters = 8
var maxIterations = 100
var seed = 25 

//Obtener colección de muestras para entrenar el K means
var trainingdataset = s2_index.sample({
  region : aoi,
  scale : 10,
  numPixels:100,
  seed : seed,
  geometries :true
});
//Segmentation
var gmeanseg = ee.Algorithms.Image.Segmentation.GMeans(s2_index)


//Aplicar el algoritmo K means
var kmeans = ee.Clusterer.wekaKMeans(numClusters).train(trainingdataset)


var clustered =s2_index.cluster(kmeans,"kmeans");
//var kmeansclus =s2_index.cluster(kmeanseg,"kmeanseg")
//print(kmeansclus,"kmeanseg")

//Ver en el mapa
Map.addLayer(clustered.randomVisualizer(),{},"K-means Clasificación")
Map.addLayer(gmeanseg.randomVisualizer(),{},"K-means Segmentation");

//Clasificación Supervisada

var trainingPolygons = ee.FeatureCollection("projects/bluemx/assets/3a_Rev_1000_points_CONABIO_2020")

trainingPolygons = trainingPolygons.filter(ee.Filter.neq('ObsID', 0));

// Crear un diccionario para reclasificación como ee.Dictionary
var reclassMapping = ee.Dictionary({
  1: 1, // Manglar
  2: 2, // Manglar perturbado
  3: 3, // Cuerpo de agua
  4: 5, // Agricola - Pecuaria
  5: 5, // Desarrollo antropico
  6: 5, // Otra vegetacion
  7: 4, // Otros humedales
  8: 5,  // Sin vegetacion
});

// Función para agregar la columna `Reclass`
var reclassify = function(feature) {
  var originalClass = feature.getNumber('ObsID'); // Obtener la clase original como número
  var newClass = reclassMapping.get(originalClass); // Obtener la clase reclasificada desde el diccionario
  return feature.set('Reclass', newClass); // Asignar la nueva clase al feature
};

// Aplicar la reclasificación a los polígonos de entrenamiento
var reclassifiedPolygons = trainingPolygons.map(reclassify);
//print('Polígonos reclasificados:', reclassifiedPolygons.first());
// 4. Dividir el conjunto de datos en entrenamiento y validación
var withRandom = reclassifiedPolygons.randomColumn('random', 40); // Semilla fija
var trainingSet = withRandom.filter(ee.Filter.lt('random', 0.6));
var validationSet = withRandom.filter(ee.Filter.gte('random', 0.4));

print(trainingPolygons, "Poligonos 1")
print(trainingSet, "Poligonos entrenamiento 1")
print(validationSet, "Poligonos validacion 1")

// 6. Extraer valores de la imagen para cada polígono
var trainingData = s2_index.sampleRegions({
  collection: trainingSet,
  properties: ['Reclass', 'PointID', "Obs2"], // Usar la nueva columna reclasificada
  scale: 30
});

print(trainingData, "Datos entrenamiento")

//Seleccionar bandas a clasificar
var selectedBands = [
        //'B2',
        'B2', 
        'B4',
        'B8', 
       'B11',
        'NDVI',
        'CMRI', 
        'MVI', 
        //'MDI', 
           'NDWI', 
        'elevation', 
           'MDI',
        'VV',
        'VH'
        ]

// 7. Entrenar el modelo
var classifier = ee.Classifier.smileRandomForest(100).train({
  features: trainingData,
  classProperty: 'Reclass',
  inputProperties: selectedBands
});


// 8. Clasificar los datos de validación
var validationData = s2_index.sampleRegions({
  collection: validationSet,
  properties: ['Reclass'], 
  scale: 30
});

// Clasificar y validar con la columna reclasificada
var validated = validationData.classify(classifier);
var errorMatrix = validated.errorMatrix('Reclass', 'classification');
print('Matriz de error:', errorMatrix);
print('Precisión general:', errorMatrix.accuracy());
print('Índice Kappa:', errorMatrix.kappa());
print('Precisión por clase:', errorMatrix.consumersAccuracy());

// 9. Clasificar la imagen completa
var classified = s2_index.classify(classifier).clip(aoi);

Map.addLayer(classified.randomVisualizer(),{},"Clasificación Supervisada")