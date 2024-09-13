.check_ggtheme <- function(ggtheme, call = rlang::caller_env()) {

  # Check if the provided theme is a suitable theme.
  if (inherits(ggtheme, "theme")) {

    # Check if the theme is complete.
    if (!attr(ggtheme, "complete")) {
      rlang::abort(
        paste0(
          "The plotting theme is not complete. The most likely cause is lack of a valid template, such as ",
          "theme_familiar or ggplot2::theme_light. Note that ggplot2::theme is designed to tweak existing ",
          "themes when creating a plot."
        ),
        call = call
      )
    }

  } else if (is.null(ggtheme)) {
    ggtheme <- ggplot2::theme_light(base_size = 9)

  } else {
    rlang::abort(
      paste0("The provided ggplot2 theme is not a theme. Found: ", class(ggtheme)),
      call = call
    )
  }

  return(ggtheme)
}
