import sys
from scipy.stats import fisher_exact

# Usage: python3 fisher_filter.py <check_file> <N_CTRL> <N_AD> <TargetGroup> <alpha>
file_path = sys.argv[1]
n_ctrl = int(sys.argv[2])
n_ad = int(sys.argv[3])
alpha = float(sys.argv[4])
beta = float(sys.argv[5])

with open(file_path, 'r') as f:
    for line in f:
        cols = line.strip().split()
        if len(cols) < 5: continue
        
        try:
            # Based on your logic:
            # k_ctrl = $3 + $5
            # k_ad   = $4
            k_ctrl = int(cols[2]) + int(cols[3])
            k_ad = int(cols[4])
            
            # Contingency Table:
            # [[Group_A_Variant, Group_A_NoVariant], [Group_B_Variant, Group_B_NoVariant]]
            table = [[k_ctrl, max(0, n_ctrl - k_ctrl)], 
                     [k_ad,   max(0, n_ad - k_ad)]]
            # Test enrichment: 
            # If target is CTRL, we test if prop_CTRL > prop_AD (alternative='greater')
            # If target is AD, we test if prop_AD > prop_CTRL (alternative='less' for the CTRL vs AD table)
            _, p_val = fisher_exact(table)

            if p_val < alpha or (k_ctrl > 0 and k_ad == 0) or (k_ad > 0 and k_ctrl == 0) :
                if k_ctrl/n_ctrl < beta and k_ad/n_ad < beta:
                    print(line.strip())

        except (ValueError, IndexError):
            continue
