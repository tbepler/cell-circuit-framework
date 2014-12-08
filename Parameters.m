classdef Parameters
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant = true);
        %timescale MINUTES, %E. coli volume ~ 1.1 ul
        DEFAULT_RNAP = 4600; % ~4600 in E. coli - http://www.ncbi.nlm.nih.gov/pubmed/22624875
        DEFAULT_RIBO = 55000; % ~55,000 in E.coli ~200,000 in yeast - http://www.ncbi.nlm.nih.gov/pubmed/22624875
        
        DEFAULT_K_TXN = 1050; % ~40-80 nt/s, 900 - 1200 bp per mRNA = ~17.5/s
        DEFAULT_K_TLN = 1050; %~20 aa/s, 300-400 aa per protein = ~17.5 proteins/s
        
        DEFAULT_PROT_DEG = 1.67e-3; % avg 0.1 / h, 0.03 - 0.82 /h reasonable,  proteosome: ~5 proteins/m - BIONUMBERS
        MAX_PROT_DEG = 1.37e-2;
        MIN_PROT_DEG = 5e-4;
        
        DEFAULT_RNA_DEG = 3.47e-2 % avg half life ~20 m ( 0.0347 /m), 3 (0.231 /m) - 90 (0.0077 /m) m reasonable - BIONUMBERS
        MAX_RNA_DEG = 2.31e-1;
        MIN_RNA_DEG = 7.7e-3;
        
        DEFAULT_RIBO_ON = 1.5398; % ~ 1.5398e-09 molecules^-1 minutes^-1 , 1.7e7 / Ms http://www.jbc.org/content/255/1/225.full.pdf
        DEFAULT_RIBO_OFF = 0.9 % 0.015 / s
        DEFAULT_RIBO_INIT = 2.4; % 0.02 - 0.06 /s http://www.jbc.org/content/289/41/28160/F7.expansion.html
        MAX_RIBO_INIT = 3.6;
        MIN_RIBO_INIT = 1.2;
        
        DEFAULT_RNAP_ON = 9.0575e-1; % strong: > 9e6 /Ms, weak: 9.6e5 /Ms http://www.ncbi.nlm.nih.gov/pmc/articles/PMC350123/pdf/pnas00497-0095.pdf
        MAX_RNAP_ON = 9.0575; %1e8 /Ms
        MIN_RNAP_ON = 8.6952e-3;
        
        DEFAULT_RNAP_OFF = 4.8e-3; %strong: 1.7e-4 /s, weak: 3.3e-5 /s http://www.ncbi.nlm.nih.gov/pmc/articles/PMC350123/pdf/pnas00497-0095.pdf
        MAX_RNAP_OFF = 5e-2;
        MIN_RNAP_OFF = 2e-4;
        
        DEFAULT_RNAP_INIT = 2; % initiation rate strong: ~0.04 /s , weak: 0.024 /s http://www.ncbi.nlm.nih.gov/pmc/articles/PMC350123/pdf/pnas00497-0095.pdf
        MAX_RNAP_INIT = 3.6;
        MIN_RNAP_INIT = 1.2;
        
        DEFAULT_TF_ON = 5e-1; % 1e7 ... 1e6 /Ms http://www.pnas.org/content/109/41/16540.full.pdf
        MAX_TF_ON = 9.0575e-1
        MIN_TF_ON = 9.0575e-2;
       
        DEFAULT_TF_OFF = 6;% 1e-2 ... 1e1 /s http://www.pnas.org/content/109/41/16540.full.pdf
        MAX_TF_OFF = 6e2;
        MIN_TF_OFF = 6e-1;
        
        DEFAULT_TF_COOP = 2; % avg 2, 1 - 10 reasonable -BIONUMBERS
        MIN_TF_COOP = 1;
        MAX_TF_COOP = 10;
        
        
        %timescale minutes
        RELAXATION_OSCILLATOR_DEFAULT_PARAMETERS = ...
        { % ~2000 proteins per mRNA in the cell
          Parameters.DEFAULT_RNAP - 3100 % RNAP ... 1
          Parameters.DEFAULT_RIBO - 40000 % Ribo ... 2
          Parameters.DEFAULT_K_TXN % k_txn ... 3 CURRENTLY UNUSED BY SYS
          Parameters.DEFAULT_K_TLN  % k_tln ... 4 CURRENTLY UNUSED BY SYS
          Parameters.DEFAULT_PROT_DEG  % x_deg ... 5
          0  % x_init ...6
          Parameters.MAX_PROT_DEG  % y_deg ...7
          0  % y_init ...8
          Parameters.DEFAULT_RNA_DEG  % x_rna_deg ...9
          Parameters.DEFAULT_RIBO_ON % x_ribo_on ...10
          Parameters.DEFAULT_RIBO_OFF  % x_ribo_off ... 11
          Parameters.DEFAULT_RIBO_INIT  % x_ribo_init ...12
          0  % x_rna_init ...13
          Parameters.DEFAULT_RNA_DEG % y_rna_deg ...14
          Parameters.DEFAULT_RIBO_ON  % y_ribo_on ...15
          Parameters.DEFAULT_RIBO_OFF  % y_ribo_off ...16
          Parameters.MAX_RIBO_INIT  % y_ribo_init ...17
          0  % y_rna_init ...18
          Parameters.MAX_RNAP_ON  % x_rnap_on ... 19
          Parameters.MIN_RNAP_OFF  % x_rnap_off ... 20
          Parameters.MAX_RNAP_INIT  % x_rnap_init ...21
          Parameters.MIN_TF_ON  % x_bindY ...22
          Parameters.MAX_TF_OFF  % x_unbindY ...23
          1  % x_copy_number ... 24
          Parameters.MIN_RNAP_ON % y_basal_rnap_on ...25
          Parameters.MAX_RNAP_OFF % y_basal_rnap_off ...26
          Parameters.MIN_RNAP_INIT % y_basal_rnap_init ...27
          Parameters.MAX_RNAP_ON  % y_bound_rnap_on ...28
          Parameters.MIN_RNAP_OFF  % y_bound_rnap_off ...29
          Parameters.MAX_RNAP_INIT  % y_bound_rnap_init ...30
          Parameters.DEFAULT_TF_ON  % y_bindX ...31
          Parameters.MAX_TF_OFF  % y_unbindX ...32
          Parameters.MAX_TF_ON  % y_bindY ...33
          Parameters.MIN_TF_OFF  % y_unbindY ...34
          1  % y_copy_number ... 35
          3 % x_coop ...36   
          6  % y_coop ...37
            };
    end
    
    
end

