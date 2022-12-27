# Step 1: Read in the data
mash.data <- mashr::mash_set_data(beta.mat, sd.mat)

print('Identifying a subset of strong tests ...')
m.1by1 <- mashr::mash_1by1(
  mash.data,
  control=list(numiter.em=200),
  optmethod='mixSQP'
)
strong.subset <- mashr::get_significant_results(m.1by1)
