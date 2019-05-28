// 测试能否实现谐波分析差值处理 
var sentinel = ee.ImageCollection("COPERNICUS/S2"),
    usda = ee.ImageCollection("USDA/NASS/CDL"),
    sample_sunflower_17 = /* color: #d63000 */ee.Geometry.MultiPoint(
        [[-100.31620454373495, 44.93244952141173],
         [-99.91461437695517, 44.86352024238593],
         [-100.29249098287397, 44.732821318882415],
         [-99.56838578272914, 44.62929137471852],
         [-100.7453479468432, 44.41871592786056],
         [-101.91365989875743, 44.28218982630365],
         [-101.26060545731707, 44.55182101229691],
         [-100.05471645850298, 44.47544920837516],
         [-99.87427593503537, 44.414482839324485],
         [-99.83523528352498, 44.22714038583697],
         [-100.01578242985158, 44.68287266870666],
         [-100.3898606343418, 46.156528028100595],
         [-100.20794482084301, 46.091288981069205],
         [-102.68819541499175, 46.41754167651417],
         [-100.57595670436808, 46.83248416284782],
         [-101.24519361700845, 45.35300632133485],
         [-101.30149854864908, 45.39207784229937],
         [-102.41933490710278, 46.16723749827688],
         [-102.44941365361564, 46.30963260964439],
         [-102.43343741396387, 46.362208372218326]]);
var usa_shape = ee.FeatureCollection("users/benkelou/usa_shp");
var qh_samples = ee.FeatureCollection("users/benkelou/SOW/SAMPLE/Qinhai_crop_sample");

/********** 预处理 ***********

/* 去云处理 
/**
 * Function to mask clouds using the Sentinel-2 QA band
 * @param {ee.Image} image Sentinel-2 image
 * @return {ee.Image} cloud masked Sentinel-2 image
 */
function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000)
    .set('system:time_start', image.get('system:time_start'))
    .set('system:time_end', image.get('system:time_end'));
  
}


//  通过置信波段usda选择作物（ 通过类似去云的方法来做）

var usda_filtered = usda
  .filterDate('2017')
  .map(function(image){
    
    var criterion = 85;
    var confidence = image.select('confidence');
    var mask = confidence.lt(criterion).eq(0);

    return image.updateMask(mask).select('cropland');
    
  });


/*
// 美国的油菜太少了，前面的挑选条件太苛刻前面的挑选条件太苛刻，放宽才行   
var usda_filtered = usda
  .filterDate('2017')
  .map(function(image){

    return image.select('cropland');
    
  });
*/
  
  
//  通过频率挑选合格作物（mode为众数）
var usda_reduced = usda_filtered
  .reduce(ee.Reducer.mode());

var crops = ee.Dictionary({
  //rapeseed: 34,
  sunflower: 6
}).toImage();

var crop_select = function(usda_reduced, crops){
  
  //var mask_rapeseed = usda_reduced.eq(crops.select('rapeseed'));
  var mask_rapeseed = usda_reduced.eq(crops.select('sunflower'));

  //var temp = ee.Image(usda_reduced.updateMask(mask_rapeseed).rename('rapeseed'));
  var temp = ee.Image(usda_reduced.updateMask(mask_rapeseed).rename('sunflower'));
  
  return temp;
  
};  
  
var crop_selected = crop_select(usda_reduced, crops);


// 执行预处理，选择数据
// 可能矢量有点问题，这个选不出来 
var boundary = ee.FeatureCollection(usa_shape)
  .filter(ee.Filter.inList('STATE_NAME',['North_Dakota','South_Dakota']));

//var boundary = ee.FeatureCollection(cn_shape). 
//  filter(ee.Filter.inList('State_Name',['Qinghai']));

// Map the function over one year of data and take the median.
// Load Sentinel-2 TOA reflectance data.
var sentinel_fliter = sentinel.filterDate('2017-01-01', '2017-12-31')
  // Pre-filter to get less cloudy granules.
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(maskS2clouds)
  .filterBounds(boundary);
  
var senrinel_for_fitting = sentinel.filterDate('2017-01-01', '2017-12-31')
  //  拿这个影像集承载谐波分析拟合的结果，补全去云处理去掉的值 
  .filterBounds(boundary)
  .map(function(image){
    return image.clip(boundary);
  });
  
// 选择波段并裁剪
var sentinel_select = sentinel_fliter.map(function(image){
  
  return image
    .select(['B2', 'B3', 'B4', 'B8'])
    .clip(boundary);
  
});
  
  
// 计算三个指标
var dataset = sentinel_select.map(function(image){
  
  var sai = image
    .expression(
      'Green / (0.6131 * Blue + 0.3869 * Red)', {
        'Blue': image.select('B2'),
        'Green': image.select('B3'),
        'Red': image.select('B4')})
    .rename('fai');
   
   var blue = image.select('B2');
   var green = image.select('B3');
   var red = image.select('B4');
   var nir = image.select('B8');
        
   var symmetry = red.subtract(blue).divide(red.add(blue)).rename('symmetry');
    //image.normalizedDifference(['B4', 'B2']).rename('symmetry');
   var indicator = red.subtract(green).divide(red.add(green)).rename('indicator');
    //image.normalizedDifference(['B4', 'B3']).rename('indixator');
    var ndvi = nir.subtract(red).divide(nir.add(red)).rename('ndvi');
    //image.normalizedDifference(['B8', 'B3']).rename('ndvi');
    
  return image
    .addBands(sai)
    .addBands(symmetry)
    .addBands(indicator)
    .addBands(ndvi);
    
});

//print(dataset);

//********** 加一个谐波分析的功能看看 **********


var harmonic_analysis = function(image, target){
  
  var timeField = 'system:time_start';
  
  var filter_image = function(inputimage){
    return inputimage
    //.filterBounds(boundary)
    // Use this function to add variables for timeseries, time and a constant
    // to imagery.
    .map(function(image) {
      // Compute time in fractional years since the epoch.
      var date = ee.Date(image.get(timeField));
      var years = date.difference(ee.Date('2017-01-01'), 'year');
      // Return the image with the added bands.
      return image
        // Add a time band.
        .addBands(ee.Image(years).rename('t').float())
        // Add a constant band.
        .addBands(ee.Image.constant(1));
    });
  };
  
  // 选择roi所在的数据并添加时间波段
  var filteredimage = filter_image(image);
  //  为输出拟合数据添加时间和常数波段    
  var computeimage = filter_image(senrinel_for_fitting); 
  
  
  //print('filteredsar', filteredsar);

  //The setup for fitting the model is to first add the harmonic variables (the third and fourth terms of equation 2) to the image collection.

  // Use these independent variables in the harmonic regression.
  var harmonicIndependents = ee.List(['constant', 't', 'cos1', 'sin1', 'cos2', 'sin2', 'cos3', 'sin3'/*, 'cos4', 'sin4', 'cos5', 'sin5', 'cos6', 'sin6'*/]);

  // Add harmonic terms as new image bands.
  var add_harmonic_paremeters = function(inputimage){
    return inputimage
      .map(function(image) {
        var timeRadians1 = image.select('t').multiply(2 * Math.PI);//第一谐波
        var timeRadians2 = image.select('t').multiply(2 * 2 * Math.PI);// 第二谐波
        var timeRadians3 = image.select('t').multiply(2 * 3 * Math.PI);// 第二谐波
        //var timeRadians4 = image.select('t').multiply(2 * 4 * Math.PI);// 第二谐波
        //var timeRadians5 = image.select('t').multiply(2 * 5 * Math.PI);// 第二谐波
        //var timeRadians6 = image.select('t').multiply(2 * 6 * Math.PI);// 第二谐波
        return image
          .addBands(timeRadians1.cos().rename('cos1'))
          .addBands(timeRadians1.sin().rename('sin1'))
          .addBands(timeRadians2.cos().rename('cos2'))
          .addBands(timeRadians2.sin().rename('sin2'))
          .addBands(timeRadians3.cos().rename('cos3'))
          .addBands(timeRadians3.sin().rename('sin3'));
          /*
          .addBands(timeRadians4.cos().rename('cos4'))
          .addBands(timeRadians4.sin().rename('sin4'))
          .addBands(timeRadians5.cos().rename('cos5'))
          .addBands(timeRadians5.sin().rename('sin5'))
          .addBands(timeRadians5.cos().rename('cos6'))
          .addBands(timeRadians5.sin().rename('sin6'));
          */
    });
  };
  
  var harmonicimage = add_harmonic_paremeters(filteredimage);
  var fittingimage = add_harmonic_paremeters(computeimage);

  //Fit the model as with the linear trend, using the linearRegression() reducer:

  // The output of the regression reduction is a 4x1 array image.
  var harmonicTrend = harmonicimage
    .select(harmonicIndependents.add(target))
    .reduce(ee.Reducer.linearRegression({
      numX: harmonicIndependents.length(), 
      numY: 1}));

  //Plug the coefficients in to equation 2 in order to get a time series of fitted values:

  // Turn the array image into a multi-band image of coefficients.
  var harmonicTrendCoefficients = harmonicTrend.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([harmonicIndependents]);

  // Compute fitted values.
  var fittedHarmonic = fittingimage.map(function(image) {
    // 计算最终的拟合值 
    return image
      .addBands(
        image.select(harmonicIndependents)
        .multiply(harmonicTrendCoefficients)
        .reduce('sum')
        .rename('fitted'));
    });
    
    return fittedHarmonic;

};

// 执行谐波分析 
var fai_analysis = harmonic_analysis(dataset, 'fai');
var symmetry_analysis = harmonic_analysis(dataset, 'symmetry');
var indicator_analysis = harmonic_analysis(dataset, 'indicator');
var ndvi_analysis = harmonic_analysis(dataset, 'ndvi');
print(fai_analysis)

//print(fai_analysis);
//print(symmetry_analysis);
//print(indicator_analysis);
//print(ndvi_analysis);

//********** 显示 **********

Map.centerObject(boundary);

//Map.addLayer(crop_selected.select('rapeseed'), {color: 'red'}, 'rapeseed');
//Map.addLayer(dataset.select('fai'), {}, 'fai');
//Map.addLayer(dataset.select('symmetry'), {}, 'symmetry');
//Map.addLayer(dataset.select('indicator'), {}, 'indicator');
//Map.addLayer(dataset.select('ndvi'), {}, 'ndvi');


Map.addLayer(crop_selected.select('sunflower'), {color: 'red'}, 'sunflower');
//Map.addLayer(fai_analysis, {}, 'fai');
//Map.addLayer(symmetry_analysis.select('fitted'), {}, 'symmetry');
//Map.addLayer(indicator_analysis.select('fitted'), {}, 'indicator');
//Map.addLayer(ndvi_analysis.select('fitted'), {}, 'ndvi');



var plot = function(image, roi, bands){
  // Plot the fitted model and the original data at the ROI.
  print(ui.Chart.image.series(
  // 这里记得改一下分辨率 
  image.select(bands), roi, ee.Reducer.mean(), 10)
    //.setSeriesNames(bands)
    .setOptions({
//      title: plot_title,
      lineWidth: 1,
      pointSize: 3,
  }));
};


//********** 获得样本点对应的数据 **********
/*
var points_dataset = dataset.getRegion(sample_sunflower_17, 10);

var property_name = points_dataset.get(0);
var under_process = points_dataset.remove(property_name);
//print(property_name);
//print(under_process);



var export_dataset = ee.FeatureCollection(under_process.map(function(el){
  el = ee.List(el); // cast every element of the list
  var log = ee.Number(el.get(1));
  //print(log)
  var lat = ee.Number(el.get(2));
  var geom = ee.Geometry.Point(log, lat);
  return ee.Feature(geom, {
    'longitude':ee.Number(el.get(1)), 
    'latitude':ee.Number(el.get(2)),
    'time':ee.Date(el.get(3)),
    'B2':ee.Number(el.get(4)),
    'B3':ee.Number(el.get(5)),
    'B4':ee.Number(el.get(6)),
    'B8':ee.Number(el.get(7)),
    'sai':ee.Number(el.get(8)),
    'symmetry':ee.Number(el.get(9)),
    'indicator':ee.Number(el.get(10)),
    'ndvi':ee.Number(el.get(11))
  });
}));
//print(export_dataset);


// 输出数据 
Export.table.toDrive({
  collection: export_dataset,
  description: 'timeseries_sunflower',
  fileFormat: 'CSV'
});

*/

//********** 绘图 **********

// 建立点击事件，选择roi点
Map.style().set('cursor','crosshair');//设置地图样式，鼠标样式为十字型
//var panel = ui.Panel({style:{width:'400px'}}).add(ui.Label('点击地图影像'));// 创建panel以存放charts
//ui.root.add(panel);//将panel加到地图上

// 建立地图点击事件
Map.onClick(function(coords) {
  
  //panel.clear();//清除画布panel
  var roi = ee.Geometry.Point(coords.lon, coords.lat);
  
  //plot(fai_analysis, roi, ['fitted']);
  //plot(symmetry_analysis, roi, ['symmetry','fitted']);
  plot(indicator_analysis, roi, ['fitted']);
  //plot(ndvi_analysis, roi, ['ndvi','fitted']);
  //plot(dataset, roi, ['fai']);
  //plot(dataset, roi, ['symmetry']);
  plot(dataset, roi, ['indicator']);
  //plot(dataset, roi, ['ndvi']);
  //plot(dataset, roi, ['SAI', 'ndvi', 'symmetry', 'indicator']);
  //plot(sentinel_select, roi, ['B2', 'B3', 'B4', 'B8']);

}); // 这个是onclick事件的结尾 
              
                  
                  










