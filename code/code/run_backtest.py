import numpy as np
import pandas as pd
from backtesting import MarkovRegimeBacktest
import matplotlib.pyplot as plt
import seaborn as sns
import os
from datetime import datetime

# Set style for plots
plt.style.use('default')
sns.set_theme()

def main():
    # Create output directory for results
    output_dir = 'backtest_results'
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize backtester
    backtester = MarkovRegimeBacktest(confidence_level=0.95)
    
    # Get current directory and construct path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    matlab_file = os.path.join(current_dir, 'matlab_results.mat')
    print(f"Loading MATLAB results from: {matlab_file}")
    backtester.load_matlab_results(matlab_file)
    
    # Generate comprehensive report with statistical tests
    report = backtester.generate_extended_summary_report()
    
    # Generate visualizations
    backtester.plot_regime_transitions(save_path=os.path.join(output_dir, 'regime_transitions.png'))
    backtester.plot_regime_distributions(save_path=os.path.join(output_dir, 'regime_distributions.png'))
    backtester.plot_qq_plots(save_path=os.path.join(output_dir, 'qq_plots.png'))
    
    # Print results
    print("\n=== Backtesting Results ===")
    print("\nRisk Measures by Regime:")
    for regime, metrics in report['risk_measures'].items():
        print(f"\n{regime.upper()}:")
        for metric, value in metrics.items():
            print(f"{metric}: {value:.4f}")
    
    print("\nViolation Analysis:")
    print(f"Observed violation ratio: {report['violation_ratio']:.4f}")
    print(f"Expected violation ratio: {report['expected_violations']:.4f}")
    
    print("\nRegime Proportions:")
    for regime, prop in report['regime_proportions'].items():
        print(f"{regime}: {prop:.2%}")
    
    print("\nStatistical Tests:")
    for regime, tests in report['statistical_tests'].items():
        if regime == 'regime_comparison':
            print(f"\nRegime Comparison Tests:")
        else:
            print(f"\n{regime.upper()} Tests:")
        
        for test_name, results in tests.items():
            print(f"{test_name}:")
            print(f"  Statistic: {results['statistic']:.4f}")
            print(f"  P-value: {results['p_value']:.4f}")

if __name__ == "__main__":
    main()
