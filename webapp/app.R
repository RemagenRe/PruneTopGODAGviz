library(shiny)
library(visNetwork)

nodes <- data.frame(
  id = c("GO:0008150", "GO:0009987", "GO:0008152", "GO:0006629",
         "GO:0044238", "GO:0071704"),
  label = c("biological_process", "cellular process", "metabolic process",
            "lipid metabolic process", "primary metabolic process",
            "organic substance metabolic process"),
  group = c("root", "term", "term", "term", "ancestor", "ancestor"),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = c("GO:0008150", "GO:0008150", "GO:0008152", "GO:0008152",
           "GO:0006629"),
  to = c("GO:0009987", "GO:0008152", "GO:0006629", "GO:0044238",
         "GO:0071704"),
  stringsAsFactors = FALSE
)

ui <- fluidPage(
  tags$head(tags$style(HTML(
    "body{font-family:Arial,Helvetica,sans-serif;background:#f4f5f7;}
     .well{background:#fff;border:1px solid #e2e4e8;border-radius:6px;}
     h2{font-size:22px;margin:10px 0 4px 0;}
     .note{color:#56606b;font-size:12px;line-height:1.5;}
     .pill{display:inline-block;padding:2px 8px;border-radius:10px;
           background:#CFEAF1;font-size:11px;margin-right:5px;}"
  ))),
  titlePanel("PruneTopGODAGviz browser preview"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Interactive DAG demo"),
      p(class = "note",
        "This is a browser preview of the interactive editor. Full annotation ",
        "upload, enrichment and export run in the R package."),
      actionButton("reset", "Reset layout", width = "100%"),
      hr(),
      h5("Selected term"),
      verbatimTextOutput("selected"),
      hr(),
      p(class = "note",
        "Colours: significant terms use a white-gold-rose scale; ancestor ",
        "terms stay white.")
    ),
    mainPanel(
      width = 9,
      div(class = "well",
        visNetworkOutput("graph", height = "68vh")
      ),
      p(class = "note",
        "Drag nodes to arrange the graph. Click a node to inspect it. ",
        "Install the R package for enrichment and publication export.")
    )
  )
)

server <- function(input, output, session) {
  output$graph <- renderVisNetwork({
    visNetwork(nodes, edges) |>
      visNodes(shape = "box", font = list(size = 14, face = "Arial")) |>
      visEdges(arrows = "to", color = "#56606b") |>
      visGroups(groupname = "root", color = list(
        background = "#FFFFFF", border = "#000000")) |>
      visGroups(groupname = "term", color = list(
        background = "#E9C46A", border = "#000000")) |>
      visGroups(groupname = "ancestor", color = list(
        background = "#FFFFFF", border = "#000000")) |>
      visOptions(highlightNearest = TRUE, nodesIdSelection = FALSE) |>
      visInteraction(dragNodes = TRUE, zoomView = TRUE) |>
      visPhysics(enabled = FALSE) |>
      visEvents(selectNode = "function(nodes) {
        Shiny.setInputValue('selected_node', nodes.nodes[0], {priority: 'event'});
      }")
  })

  observeEvent(input$reset, {
    visNetworkProxy("graph") |> visSetSelection(nodes = NULL)
  })

  output$selected <- renderText({
    if (is.null(input$selected_node)) return("Click a node")
    i <- match(input$selected_node, nodes$id)
    if (is.na(i)) return(input$selected_node)
    paste0(nodes$id[i], "\n", nodes$label[i])
  })
}

shinyApp(ui, server)
