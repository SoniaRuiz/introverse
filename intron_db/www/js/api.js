function goAFunction(table, intronId) {
  //console.log(table.tables().nodes().to$().attr('id'));
  console.log(decodeURIComponent(intronId))
  table.search(decodeURIComponent(intronId)).draw();
  //table.column(11).search(decodeURIComponent(intronId)).draw();
  
  //table.search(decodeURIComponent(intronId)).draw();
  //alert(intronId);
  
  $('#goB').show();
  $('#goA').hide();
}

function goBFunction(tableIntrons, tableNovel) {
  //console.log(table.tables().nodes().to$().attr('id'));
  //console.log(tableNovel.tables().nodes().to$().attr('id'));
  
  
  tableIntrons.search('').draw();
  
  tableNovel.clear().destroy();
  $("#intronGeneDetail_tab1").empty();

  
  $('#goB').hide();
  $('#goA').show();
}

