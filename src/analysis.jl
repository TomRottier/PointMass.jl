# terminal speed
terminal_speed(p) = 1 / sqrt(min_CD(p))

# steady state values
steady_state_ν(α, p) = sqrt(1 / sqrt(CD(α, p)^2 + CL(α, p)^2))
steady_state_θ(α, p) = -atan(CD(α, p), CL(α, p))

# alpha giving max CL/CD
min_glide_slope(p) = atan(sqrt(min_CD(p) / (max_CD(p))))