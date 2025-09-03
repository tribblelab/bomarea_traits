# Bomarea Trait Analysis (RevBayes, bash, R)

## Step 1: Using raw data to make a nexus file of inflorescence types for Bomarea
1. Run `bomareacode.R` on `bomarea traits.xlsx`
   - This code takes the measurements of herbarium specimens and data from the Excel sheet and categorizes traits: size, sparsity, branchiness, inflorescence type (binary and tertiary).
2. You should have `.nexus` files in /data. In this example, inflorescence type will be used

## Step 2: Generate a tree using `type.nexus` and `bom_only_MAP.tre`
1. Set working directory in Zsh
   ```bash
   cd ~/Desktop/bomarea_traits/
   ```
2. Open RevBayes in Zsh and execute the script
   - Enter `rb` in your terminal to run RevBayes.
   - Run the `infl_type_ard.Rev` script until `->`.

3. In bash, combine and remove the first 10% of data:
   - In another terminal window, change directory to `~/Desktop/bomarea_traits/output`.
   - Run the code to combine the runs (`.txt` files):
     ```bash
     awk 'FNR == 1 && NR != 1 { next } { print }' size_ard_states_run_1.txt size_ard_states_run_2.txt > size_ard_states_combined.txt
     ```

4. Return to the RevBayes window (where you ran the MCMC):
   - Execute the two lines after `->`.

## Step 3: Plot the tree in R
1. Run `plotresults.R` with `infl_type_ase_ard.tree`

## Step 4: Create a Violin Plot using rates
1. Run `plotviolinrates.R` and save the plot as a `.png`.
   - This will create a violin plot showing the rates of transitions between inflorescence types.

---
# Notes

To run the script in RevBayes:
```bash
< infl_type_ard_binary.Rev
```

This is the same as Step 2, point 3, but with different file names:
1. Run the code to combine runs (`.txt` files):
   ```bash
   awk 'FNR == 1 && NR != 1 { next } { print }' branchiness_run_1.txt branchiness_run_2.txt > branchiness_combined.txt
   ```
2. Run the code to remove the first 10% of the data:
   ```bash
   total_lines=$(wc -l < branchiness_combined.txt)
   skip_lines=$((total_lines / 10))
   awk -v skip="$skip_lines" 'NR > skip || NR == 1' branchiness_combined.txt > branchiness_combined_trimmed.txt
   ```
