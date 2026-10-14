export parse_geom_file
export load_resonance_grid

"""
    parse_geom_file(path; comment_char='#')

Parse a space-separated N-dimensional geometries file shaped like

    #geom   coord1   coord2   coord3
        1   +0.0     -0.5     +0.3
        2   +0.1     -0.5     +0.3
        ...

Returns:
- `coord_names :: Vector{String}`: column names after "geom"
- `coords :: Matrix{Float64}`: matrix of size `ngeoms x ncoords`. Row `i` is geometry `i`. Column `j` is `coord_names[j]`

Enforces geometry contiguitity: sorts the coordinates by the index if the geometries are non-contiguous
"""
function parse_geom_file(
      path :: AbstractString
    ; comment_char :: Char = '#'
  )

  isfile(path) || error("Geometry file not found: $path")

  coord_names = String[]
  rows = Tuple{Int, Vector{Float64}}[]
  header_found = false
  ncols = 0

  open(path) do io

    for (linenumber, raw) in enumerate(eachline(io))

      line = strip(raw)
      isempty(line) && continue

      if startswith(line, comment_char)
        stripped = strip(lstrip(line, comment_char))
        tokens = split(stripped)

        if !header_found && length(tokens) >= 2 && lowercase(tokens[1]) == "geom"
          coord_names = String.(tokens[2:end])
          ncols = length(coord_names)
          header_found = true
        end

        continue

      end

      tokens = split(line)
      header_found || error("No header line found before data in $path. "*
                           "I need to know what to call the columns.")
      length(tokens) == 1 + ncols ||
        error("Geometry file $path, line $linenumber: expected $(1+ncols) columns, got $(length(tokens)) tokens")

      igeom = try
        parse(Int, tokens[1])
      catch
        error("Geometry file $path, line $linenumber: can't parse geometry index $(tokens[1])")
      end

      vals = Vector{Float64}(undef, ncols)

      for i in 1:ncols
        vals[i] = try
          parse(Float64, tokens[1 + i])
        catch
          error("Geometry file $path, line $linenumber: can't parse $(coord_names[i]) value $(tokens[1+i])")
        end
      end

      push!(rows, (igeom, vals))
    end
  end

  header_found  || error("Never found a recognizable header in $path. ")
  isempty(rows) && error("Geometry file $path has a header, but no data rows")

  # -- sort by geometry index
  sort!(rows; by = r -> r[1])

  # -- check for dupes and contiguity
  rows[1][1] == 1 || error("Geometry 1 is not 1 in file $path after sort: $(rows[1][1])")
  for k in 2:length(rows)
    rows[k][1] == rows[k-1][1] && error("Geometry file $path: duplicate geometry index $(rows[k][1])")
    rows[k][1] == rows[k-1][1] + 1 && continue
    error("Geometry file $path: non-contiguous elements:
          $(rows[k][:])
          $(rows[k-1][:])"
    )
  end

  # -- pack up the outputs
  ngeoms = length(rows)
  coords = Matrix{Float64}(undef, ngeoms, ncols)
  for (k, (idx, vals)) in enumerate(rows)
    @inbounds coords[k, :] .= vals
  end

  return coord_names, coords

end

"""
    load_resonance_grid(root; geomfile, fname_regex, recursive, comments, unitsE, res_detect)
      → ResonanceGrid

Loads an N-D resonance scan from a UKRmol+output directory.

Arguments:
- `root`: directory path containing geometries file and eigenphase outputs
- `geomfile`: path to geometries file
- `fname_regex`: regex for eigenphase filenames
- `res_detect`: NamedTuple of kwargs to get forwarded to `find_resonances`
- `recursive`, `comments`, `unitsE`: kwargs forwarded to `load_eigenph_records`

Returns the assembled `ResonanceGrid`
"""
function load_resonance_grid(
      root        :: AbstractString
    ; geomfile    :: AbstractString = "geometries"
    , fname_regex :: Regex = r"eigenph\.all\.geom\d+$"
    , recursive   :: Bool = true
    , comments    :: Bool = true
    , unitsE      :: Symbol = :eV
    , res_detect  :: NamedTuple = (cols = :all , Emin = :auto , min_prominence = 0.02 , smoothwin = 1)
  )

  # -- parse N-D geometries file
  gf = geomfile === :auto ? joinpath(root, "geometries") : geomfile
  isfile(gf) || error("Geometries file not found: $gf")
  coord_names, coords = parse_geom_file(gf)

  # -- load eigenphase records via
  records = load_eigenph_records(
      root
    ; geomfile = nothing
    , attach_Q = false
    , attach_geom = true
    , recursive = recursive
    , fname_regex = fname_regex
    , comments = comments
    , unitsE = unitsE
  )

  ngeoms = size(coords, 1)
  length(records) == ngeoms || error("Number of geoms ($ngeoms) does not match number of records ($(length(records)))")

  # -- find res per rec
  idxres = find_resonances_across_records(records; res_detect...)
  resonances = resonances_by_record(idxres, length(records))

  eigenph = [r.data for r in records]
  paths   = [r.path for r in records]

  # -- pack into ResonanceGrid
  return ResonanceGrid(
       coord_names
    , coords
    , eigenph
    , resonances
    , paths
    , String(root)
  )

end
