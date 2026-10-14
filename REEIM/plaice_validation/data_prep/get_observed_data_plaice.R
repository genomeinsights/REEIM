## Contemporary observed data for the North Sea plaice hold-out validation of REEIM.
##
## Builds the same objects as empirical_data/R/get_observed.data.R does for perch:
##   age_structure_obs   data.table(age, size)   individual fish, size in mm
##   size_structure_obs  data.table(L)           individual fish, 1-cm length class (floor, cm)
## plus a maturity ogive used to set LT_obs.
##
## Sources (ICES DATRAS, downloaded with icesDatras):
##   * BTS (Beam Trawl Survey), quarter 3, 2015-2020: length frequencies (HL) and
##     age-length samples (CA). Beam trawl with a fine-mesh liner, so it catches juveniles,
##     playing the role the Swedish gillnet monitoring plays for perch.
##   * NS-IBTS, quarter 1, 2015-2020: maturity (CA). Plaice spawn in Q1, so Q3 maturity
##     staging is unreliable and Q1 is used for the ogive.
## Restricted to the North Sea: shoot latitude > 51.0 N (drops the eastern Channel, 7.d) and not
## west of Great Britain (the UK BTS also fishes the Irish Sea and Bristol Channel, a different
## plaice stock: hauls south of 57 N and west of 2.5 W are dropped).
## Beam widths differ between vessels (4, 7 and 8 m), so catches are standardised per metre of
## beam per hour, a swept-area proxy (towing speeds are all close to 4 knots).
##
## Run from the repository root (needs internet access to datras.ices.dk):
##   Rscript REEIM/plaice_validation/data_prep/get_observed_data_plaice.R
##
## Writes REEIM/plaice_validation/data/observed_data_plaice_female.RData (used by
## REEIM/input/plaice/abc_settings.csv), observed_data_plaice_pooled.RData and summary CSVs.
## Raw DATRAS downloads are cached there as raw_*.rds; with the cache present no internet is needed.

suppressPackageStartupMessages({
  library(data.table)
})

years    <- 2015:2020
aphia    <- 127143      # Pleuronectes platessa (WoRMS AphiaID)
out_dir  <- "./REEIM/plaice_validation/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## DATRAS length classes: LngtCode "." and "0" are in mm, "1" and "5" in cm
to_mm <- function(cls, code) ifelse(code %in% c("1", "5"), cls * 10, cls)

## DATRAS column names differ between download formats (e.g. AgeRings vs Age).
## Rename the alternatives to the names used below, and stop with the actual names if one is missing.
alt_names <- list(
  AgeRings   = c("AgeRings", "Age", "IndividualAge"),
  CANoAtLngt = c("CANoAtLngt", "NoAtALK", "NumberAtLength"),
  HLNoAtLngt = c("HLNoAtLngt", "NumberAtLength"),
  LngtClass  = c("LngtClass", "LengthClass"),
  LngtCode   = c("LngtCode", "LengthCode"),
  SubFactor  = c("SubFactor", "SubsamplingFactor"),
  HaulDur    = c("HaulDur", "HaulDuration"),
  HaulVal    = c("HaulVal", "HaulValidity"),
  Maturity   = c("Maturity", "IndividualMaturity"),
  Sex        = c("Sex", "IndividualSex"))
need <- list(HH = c("HaulVal", "HaulDur", "ShootLat", "DataType"),
             HL = c("HLNoAtLngt", "SubFactor", "LngtClass", "LngtCode"),
             CA = c("AgeRings", "CANoAtLngt", "LngtClass", "LngtCode", "Sex", "Maturity"))
normalise <- function(d, record) {
  for (target in names(alt_names)) {
    if (target %in% names(d)) next
    hit <- intersect(alt_names[[target]], names(d))[1]
    if (!is.na(hit)) setnames(d, hit, target)
  }
  miss <- setdiff(need[[record]], names(d))
  if (length(miss))
    stop(record, " record lacks ", paste(miss, collapse = ", "),
         ". Columns present: ", paste(names(d), collapse = ", "))
  d
}

get <- function(record, survey, quarter) {
  f <- file.path(out_dir, sprintf("raw_%s_%s_Q%d_%d-%d.rds", record, survey, quarter, min(years), max(years)))
  if (file.exists(f)) d <- readRDS(f)
  else {
    if (!requireNamespace("icesDatras", quietly = TRUE)) stop("install.packages(\"icesDatras\") to download ", f)
    d <- as.data.table(icesDatras::getDATRAS(record = record, survey = survey, years = years, quarters = quarter))
    saveRDS(d, f)
  }
  normalise(as.data.table(d), record)
}

species_filter <- function(d) {
  if ("Valid_Aphia" %in% names(d)) d[Valid_Aphia == aphia]
  else d[SpecCodeType == "W" & SpecCode == aphia]
}

haul_key <- c("Survey", "Quarter", "Country", "Ship", "Gear", "StNo", "HaulNo", "Year")

## ---- Hauls -----------------------------------------------------------------
hh <- get("HH", "BTS", 3)
hh <- hh[HaulVal == "V" & ShootLat > 51.0 & !(ShootLat < 57 & ShootLong < -2.5)]
hh[, beam_m := as.numeric(sub("^BT([0-9]+).*$", "\\1", Gear))]
if (anyNA(hh$beam_m)) stop("Unknown beam width for gears: ", paste(unique(hh[is.na(beam_m), Gear]), collapse = ", "))
print(hh[, .N, by = .(Country, Gear, beam_m)])
hh <- unique(hh[, c(haul_key, "HaulDur", "ShootLat", "ShootLong", "StatRec", "DataType", "beam_m"), with = FALSE])

## ---- Length distribution (HL): numbers per hour, summed over hauls ----------
hl <- species_filter(get("HL", "BTS", 3))
hl <- merge(hl, hh, by = haul_key, suffixes = c("", ".hh"))
hl <- hl[!is.na(LngtClass) & HLNoAtLngt > 0]
hl[, L_mm := to_mm(LngtClass, LngtCode)]
## DataType R: raw counts, raise by SubFactor and standardise to one hour.
## DataType C: already numbers per hour. Then divide by beam width (m).
hl <- hl[SubFactor > 0]   # -9 = missing
hl[, cpue := fifelse(DataType == "C", HLNoAtLngt, HLNoAtLngt * SubFactor * 60 / HaulDur) / beam_m]
hl[, L_cm := floor(L_mm / 10)]
len_freq <- hl[, .(cpue = sum(cpue), n_measured = sum(HLNoAtLngt)), by = .(Year, L_cm)]
fwrite(len_freq[order(Year, L_cm)], file.path(out_dir, "BTS_Q3_length_frequency.csv"))

## ---- Age-length samples (CA) ---------------------------------------------
ca <- species_filter(get("CA", "BTS", 3))
ca <- merge(ca, hh[, haul_key, with = FALSE], by = haul_key)
ca <- ca[!is.na(AgeRings) & AgeRings >= 0 & !is.na(LngtClass)]
ca[, L_mm := to_mm(LngtClass, LngtCode)]
ca[, L_cm := floor(L_mm / 10)]
ca[, n := fifelse(is.na(CANoAtLngt) | CANoAtLngt < 1, 1, CANoAtLngt)]

## Raw-otolith version by sex, for checking dimorphism (biased by length stratification)
fwrite(ca[, .(n = sum(n), mean_L_mm = weighted.mean(L_mm, n)), by = .(Sex, age = AgeRings)][order(Sex, age)],
       file.path(out_dir, "BTS_Q3_raw_mean_length_at_age_by_sex.csv"))

## ---- Female share at length --------------------------------------------------
## BTS length frequencies are not sexed, so the female length distribution is the CPUE at
## length times the female share at that length in the otolith samples (years pooled).
## Classes with fewer than 20 sexed fish borrow the share of the nearest class that has 20.
sexed  <- ca[Sex %in% c("F", "M"), .(nF = sum(n[Sex == "F"]), nS = sum(n)), by = L_cm][order(L_cm)]
ok     <- sexed[nS >= 20]
pF_at  <- function(l) ok$nF[which.min(abs(ok$L_cm - l))] / ok$nS[which.min(abs(ok$L_cm - l))]
len_freq[, pF := sapply(L_cm, pF_at)]
fwrite(unique(len_freq[, .(L_cm, pF)])[order(L_cm)], file.path(out_dir, "BTS_Q3_female_share_at_length.csv"))

## ---- Observed-data sets --------------------------------------------------------
## "female": female length distribution and female age-length key (default for REEIM, since
##           L_T is compared with female PMRNs and LT_obs is the female L50).
## "pooled": both sexes (and unsexed fish), as was done for perch.
## BTS otoliths are length-stratified, so raw mean size-at-age is biased. The age-length key
## (per year) is applied to the CPUE length distribution instead.
build_set <- function(set) {
  lf  <- copy(len_freq)
  cas <- if (set == "female") ca[Sex == "F"] else ca
  if (set == "female") { lf[, cpue := cpue * pF]; lf[, n_measured := n_measured * pF] }

  ## Individual-fish table for the KDE comparison. density() in run_REEIM_ABC.R is unweighted,
  ## so CPUE at length is turned into integer counts, keeping the total near the number measured.
  lf_all <- lf[, .(cpue = sum(cpue)), by = L_cm]
  scale  <- sum(lf$n_measured) / sum(lf_all$cpue)
  size_structure_obs <- lf_all[, .(L = rep(L_cm, round(cpue * scale)))]

  alk <- cas[, .(n = sum(n)), by = .(Year, L_cm, age = AgeRings)]
  alk[, p := n / sum(n), by = .(Year, L_cm)]
  naa <- merge(lf[, .(Year, L_cm, cpue)], alk[, .(Year, L_cm, age, p)],
               by = c("Year", "L_cm"), allow.cartesian = TRUE)
  ## Length classes with no otoliths in that year are dropped (share reported below).
  unaged <- lf[!alk, on = .(Year, L_cm), sum(cpue)] / sum(lf$cpue)
  naa_cm <- naa[, .(num = sum(cpue * p)), by = .(age, L_cm)]
  fwrite(naa_cm[order(age, L_cm)], file.path(out_dir, sprintf("BTS_Q3_numbers_at_age_length_%s.csv", set)))

  ## Individual (age, size in mm) pairs at the class midpoint, as for perch (mm).
  n_target <- 20000
  naa_cm[, k := round(num / sum(num) * n_target)]
  age_structure_obs <- naa_cm[k > 0, .(age = rep(age, k), size = rep(L_cm * 10 + 5, k))]

  cat(sprintf("[%s] fish: %.0f  otoliths: %.0f  CPUE share in unaged length classes: %.4f  rows size/age: %.0f / %.0f\n",
              set, sum(lf$n_measured), sum(cas$n), unaged, nrow(size_structure_obs), nrow(age_structure_obs)))
  list(age_structure_obs = age_structure_obs, size_structure_obs = size_structure_obs)
}
sets <- lapply(c(female = "female", pooled = "pooled"), build_set)


## ---- Maturity ogive (NS-IBTS Q1) --------------------------------------------
hh1 <- get("HH", "NS-IBTS", 1)
## The stock is the North Sea plus the Skagerrak; drop Kattegat hauls (east of 10.4 E, south of 57.75 N).
hh1 <- unique(hh1[HaulVal == "V" & !(ShootLong > 10.4 & ShootLat < 57.75), haul_key, with = FALSE])
mat <- species_filter(get("CA", "NS-IBTS", 1))
mat <- merge(mat, hh1, by = haul_key)
mat[, L_mm := to_mm(LngtClass, LngtCode)]
mat[, n := fifelse(is.na(CANoAtLngt) | CANoAtLngt < 1, 1, CANoAtLngt)]
## Maturity codes: 1 / 61 = immature; 2-4 and 62-66 = mature (maturing, spawning, spent, resting).
## Check the codes present before trusting this.
print(mat[, .N, by = Maturity][order(Maturity)])
mat[, mature := fcase(Maturity %in% c("1", "61", "A"), 0L,
                      Maturity %in% c("2", "3", "4", "5", "62", "63", "64", "65", "66", "B", "C", "D", "E"), 1L,
                      default = NA_integer_)]
mat <- mat[!is.na(mature) & Sex %in% c("F", "M")]
ogive <- rbindlist(lapply(c("F", "M"), function(s) {
  m <- glm(mature ~ L_mm, family = binomial, weights = n, data = mat[Sex == s])
  data.table(Sex = s, L50_mm = -coef(m)[1] / coef(m)[2],
             L25_mm = (qlogis(0.25) - coef(m)[1]) / coef(m)[2],
             L75_mm = (qlogis(0.75) - coef(m)[1]) / coef(m)[2],
             n = sum(mat[Sex == s, n]))
}))
fwrite(ogive, file.path(out_dir, "IBTS_Q1_maturity_ogive.csv"))

fwrite(mat[, .N, by = .(Sex, Maturity)][order(Sex, Maturity)], file.path(out_dir, "IBTS_Q1_maturity_codes.csv"))

for (set in names(sets)) {
  age_structure_obs  <- sets[[set]]$age_structure_obs
  size_structure_obs <- sets[[set]]$size_structure_obs
  save(age_structure_obs, size_structure_obs, ogive,
       file = file.path(out_dir, sprintf("observed_data_plaice_%s.RData", set)))
}

cat("\nHauls:", nrow(hh), "\n")
print(ogive)