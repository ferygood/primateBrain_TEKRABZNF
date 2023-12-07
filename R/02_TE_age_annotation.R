# TE annotation file
# Read the file
file_path <- "~/Downloads/Dfam_curatedonly.embl.gz"
file_lines <- readLines(gzfile(file_path))

# Initialize empty vectors to store ID, NM, and OS
ID <- character()
NM <- character()
OS <- character()

# Initialize variables to track the current chunk
current_ID <- NULL
current_NM <- NULL
current_OS <- NULL

# Define the pattern to match in the OS string
os_pattern <- "Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Homo sapiens"

# Iterate through the lines and extract information from each chunk
for (i in seq_along(file_lines)) {
    line <- file_lines[i]

    if (grepl("^ID\\s+", line)) {
        # Start of a new chunk, save the previous chunk information
        if (!is.null(current_ID) && !is.null(current_OS)) {
            ID <- c(ID, current_ID)
            NM <- c(NM, current_NM)
            OS <- c(OS, current_OS)
        }

        # Reset the variables for the new chunk
        current_ID <- gsub("^ID\\s+", "", line)
        current_NM <- NULL
        current_OS <- NULL
    } else if (grepl("^NM\\s+", line)) {
        current_NM <- gsub("^NM\\s+", "", line)
    } else if (grepl("^OS\\s+", line)) {
        os_content <- gsub("^OS\\s+", "", line)

        # Check if the OS content matches the pattern
        if (grepl(os_pattern, os_content, ignore.case = TRUE)) {
            current_OS <- os_content
        }
    }
}

# Save the last chunk information
if (!is.null(current_ID) && !is.null(current_OS)) {
    ID <- c(ID, current_ID)
    NM <- c(NM, current_NM)
    OS <- c(OS, current_OS)
}

# Create a dataframe with ID, NM, and OS
result_df <- data.frame(ID = ID, NM = NM, OS = OS, stringsAsFactors = FALSE)
write.csv(result_df, file = "~/github/pBrain/data/Dfam_TE_simiiformes.csv")
