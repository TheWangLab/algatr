ascii_alligator <- function(i) {
  cat(crayon::green(alligator(i)), sep = "\n")
}

alligator <- function(x) {
  a <- c(
    "",
    "                                 .lxk:    ",
    "                                'x0kXXc   ",
    "                             .;dXMMMMMO.  ",
    "                         .,lkXWMMMWNWO'   ",
    paste0("   .;clool:.      ..,clox0NMMMWXXXo", crayon::white("\\/")),
    paste0("  dXWK", crayon::black("xx"), "WKKKxoodkOKNWMMMMMMXO0k", crayon::white("\\/")),
    paste0("  MM0", crayon::black("xx"), crayon::white("O"), crayon::black("x"), "MMMMMMMMMMMMMWOdkc", crayon::white("\\/")),
    paste0(
      "  MMX", crayon::black("xxxx"), "MMMMMMMWMXdldkc", crayon::white("\\/"), "      ",
      crayon::bold(crayon::red(paste0(x, " done!")), "         ")
    ),
    paste0("  MMMMMMMMMMMN0xooo;", crayon::white("\\/")),
    paste0("  MMMMMMMMMMWOoc", crayon::white("/\\ /\\ /\\ /\\ /\\ /\\ /\\ /\\ /\\")),
    "  MMAPBMMMMMMMMMMEACMMMWWMWWWMIJWWWWWWKXWX",
    paste0(crayon::blue("~~~ "), "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW ", crayon::blue("~~~")),
    crayon::blue(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
    crayon::blue("  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  "),
    ""
  )
  return(a)
}
