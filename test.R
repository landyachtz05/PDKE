# would normally do a print or cat w/o using stderr, but it's the only thing
# that gets the shiny server on prod to actually print something.
cat(file=stderr(), "In test", "\n")

# Get the command-line arguments passed from the main script
args <- commandArgs(trailingOnly = TRUE)
cat(file=stderr(), "args:", length(args), "\n")

for (a in args) {
  cat(file=stderr(), "arg_v2: ", a, "\n")
  cat(file=stderr(), "arg_v2: ", a, "\n")
}
