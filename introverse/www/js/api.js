function goAFunction(table, intronId) {
  

  //console.log(table.tables().nodes().to$().attr('id'));
  //console.log(decodeURIComponent(intronId))
  table.search(decodeURIComponent(intronId)).draw();
  //table.column(11).search(decodeURIComponent(intronId)).draw();
  
  //table.search(decodeURIComponent(intronId)).draw();
  //alert(intronId);
  
  $('#goB').show();
  $('#goA').hide();
}

function goBFunction(tableIntrons, tableNovel) {
  
  tableIntrons.search('').draw();
  
  tableNovel.clear().destroy();
  $("#intronGeneDetail_tab1").empty();
  
  $('#goB').hide();
  $('#goA').show();
}

/*******************************************/
/** FUNCTIONS TO CONTROL THE BS-MODAL IMG **/
/*******************************************/

$(function(){
 $(this).on("show.bs.modal", function(e) {
   //alert("Hi")
   $("#modalVisualiseTranscriptNovel_tab1").find("div.load-container").removeClass("shiny-spinner-hidden");
   $("#modalVisualiseTranscriptNovel_tab1").find("div.load-container").addClass("shiny-spinner");

 });
 $(this).on("shown.bs.modal", function(e) {
   //alert("Hi")
   $("#modalVisualiseTranscriptNovel_tab1").find("div.load-container").addClass("shiny-spinner-hidden");
   $("#modalVisualiseTranscriptNovel_tab1").find("div.load-container").removeClass("shiny-spinner");

 });
 $(this).on("hidden.bs.modal", function(e) {
   $("#modalVisualiseTranscriptNovel_tab1 > img").remove();
 });
})



/**********************************************************/
/**** REDIRECT PARENT (in order to avoid iframe issue) ****/
/**********************************************************/

$( document ).ready(function() {
  
  var path = $("#shinyframe", window.parent.document).attr("src");
  if(path != undefined)
    window.top.location.href = path;
    
});

$(function(){
  var link = parent.document.createElement('link');
  link.type = 'image/ico';
  link.rel = 'shortcut icon';
  link.href = 'https://raw.githubusercontent.com/SoniaRuiz/introverse/main/introverse/www/introverse.ico';
  parent.document.getElementsByTagName('head')[0].appendChild(link);
  
})