# Utility script to format species names from Metaphlan profiles (For eg. s__Prevotella_melaninogenica) into italics(Prevotella~melaninogenica)
library(stringr)
library(tidyverse)

format_species_label <- function(x) {
  x %>%
    str_replace_all("^s__", "") %>%                    
    str_replace_all("_", " ") %>%                      
    str_replace("^ ", "~") %>%                   # Replace first space with ~
#    str_replace("~(sp|bacterium|unclassified)", "~\"\\1") %>%  # Start quoting if species is a multi-word like sp or bacterium
#    str_replace("(sp|bacterium|unclassified).*$", "\\0\"") %>% # Close quote at end
    paste0("italic(\"", ., "\")")                                   # Wrap in italic()
}
