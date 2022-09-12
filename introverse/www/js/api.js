function goAFunction(tableIntrons, intronId) {
  // Collapse all rows before displaying novel junctions
 /* tableIntrons.rows().every(function(){
      // If row has details expanded
      if(this.child.isShown()){
          // Collapse row details
          this.child.hide();
          $(this.node()).removeClass('shown');
      }
  });*/
  tableIntrons.search(decodeURIComponent(intronId)).draw();

  
  $('#goB').show();
  $('#goA').hide();
}

function goBFunction(tableIntrons, tableNovel) {
  /*tableIntrons = table.closest("table").DataTable()
  tableNovel = table.closest("div[data-value='Splicing Activity']").find("#intronGeneDetail_tab1").find("table").DataTable()*/
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
  //$("#modalVisualiseTranscriptNovel_tab1 > img").removeClass("d-none");
  setTimeout(function () {
      $("#modalVisualiseTranscriptNovel_tab1 > img").removeClass("d-none");
  }, 500);

 });
 $(this).on("hidden.bs.modal", function(e) {
   $("#modalVisualiseTranscriptNovel_tab1 > img").addClass("d-none");
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