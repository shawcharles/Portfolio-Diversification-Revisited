import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from scipy.io import loadmat
from typing import Dict, Tuple, List
import seaborn as sns
from datetime import datetime

class MarkovRegimeBacktest:
    def __init__(self, confidence_level: float = 0.95):
        """
        Initialize backtesting framework for Markov-modulated Lévy process model
        
        Parameters:
        -----------
        confidence_level : float
            Confidence level for VaR and CVaR calculations
        """
        self.confidence_level = confidence_level
        self.returns = None
        self.regimes = None
        
    def load_matlab_results(self, matlab_file: str) -> None:
        """
        Load results from MATLAB analysis
        
        Parameters:
        -----------
        matlab_file : str
            Path to MATLAB .mat file containing results
        """
        try:
            data = loadmat(matlab_file)
            
            # Load returns data
            self.returns = data['Data'].flatten()
            
            # Combine regime data from DataState1 and DataState2
            state1_data = data['DataState1'].flatten()
            state2_data = data['DataState2'].flatten()
            
            # Create regime indicators array
            self.regimes = np.ones(len(self.returns), dtype=int)  # Default to regime 1
            
            # Find indices where data exists in each state
            state1_valid = state1_data != 0
            state2_valid = state2_data != 0
            
            # Map the data to time points
            state1_returns = state1_data[state1_valid]
            state2_returns = state2_data[state2_valid]
            
            # Create masks for each regime
            regime1_mask = np.isin(self.returns, state1_returns)
            regime2_mask = np.isin(self.returns, state2_returns)
            
            # Assign regimes
            self.regimes[regime2_mask] = 2
            
            print(f"Loaded {len(self.returns)} observations")
            print(f"Regime 1 (Base): {np.sum(self.regimes == 1)} observations")
            print(f"Regime 2 (Spike): {np.sum(self.regimes == 2)} observations")
            
            # Verify regime proportions match P_Rt_1 and P_Rt_2
            p_rt_1 = data['P_Rt_1'][0,0]
            p_rt_2 = data['P_Rt_2'][0,0]
            actual_p1 = np.mean(self.regimes == 1)
            actual_p2 = np.mean(self.regimes == 2)
            
            print(f"\nRegime Probabilities:")
            print(f"Expected P(Regime 1): {p_rt_1:.4f}, Actual: {actual_p1:.4f}")
            print(f"Expected P(Regime 2): {p_rt_2:.4f}, Actual: {actual_p2:.4f}")
            
        except Exception as e:
            print(f"Error loading MATLAB data: {e}")
            raise
            
    def calculate_traditional_var(self, returns: np.ndarray) -> float:
        """Calculate traditional Value at Risk"""
        return np.percentile(returns, (1-self.confidence_level)*100)

    def calculate_traditional_cvar(self, returns: np.ndarray) -> float:
        """Calculate traditional Conditional Value at Risk"""
        var = self.calculate_traditional_var(returns)
        return returns[returns <= var].mean()
    
    def calculate_regime_risk_measures(self) -> Dict:
        """Calculate risk measures for each regime"""
        risk_metrics = {}
        for regime in [1, 2]:  # Base and Spike regimes
            regime_returns = self.returns[self.regimes == regime]
            risk_metrics[f'regime_{regime}'] = {
                'var': self.calculate_traditional_var(regime_returns),
                'cvar': self.calculate_traditional_cvar(regime_returns),
                'volatility': np.std(regime_returns),
                'skewness': stats.skew(regime_returns),
                'kurtosis': stats.kurtosis(regime_returns),
                'mean': np.mean(regime_returns),
                'sample_size': len(regime_returns)
            }
        return risk_metrics
    
    def calculate_violation_ratios(self, window_size: int = 252) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate VaR violation ratios using rolling windows
        
        Parameters:
        -----------
        window_size : int
            Size of rolling window in days
            
        Returns:
        --------
        tuple : (violation_ratios, expected_violations)
        """
        violations = np.zeros_like(self.returns)
        expected_violations = (1 - self.confidence_level) * np.ones_like(self.returns)
        
        for i in range(window_size, len(self.returns)):
            window_returns = self.returns[i-window_size:i]
            var = self.calculate_traditional_var(window_returns)
            violations[i] = self.returns[i] < var
            
        return violations, expected_violations
    
    def kupiec_test(self, violations: np.ndarray) -> Dict:
        """
        Perform Kupiec test for VaR violations
        
        Returns:
        --------
        dict : Test statistics and p-value
        """
        T = len(violations)
        N = np.sum(violations)
        p = 1 - self.confidence_level
        
        if N == 0:
            return {'statistic': np.nan, 'p_value': np.nan}
            
        likelihood_ratio = -2 * (np.log(((1-p)**T) * (p**N)) - 
                               np.log(((1-N/T)**T) * ((N/T)**N)))
        p_value = 1 - stats.chi2.cdf(likelihood_ratio, df=1)
        
        return {'statistic': likelihood_ratio, 'p_value': p_value}
    
    def plot_regime_transitions(self, save_path: str = None) -> None:
        """Plot regime transitions and corresponding returns"""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
        
        # Plot returns
        ax1.plot(self.returns, alpha=0.7, label='Returns', color='blue')
        ax1.set_title('Asset Returns')
        ax1.set_ylabel('Returns')
        ax1.legend()
        
        # Plot regime indicators
        ax2.plot(self.regimes, label='Regime', color='red', linewidth=1)
        ax2.set_title('Regime Indicators (1: Base, 2: Spike)')
        ax2.set_ylabel('Regime')
        ax2.set_xlabel('Time')
        ax2.legend()
        
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        plt.show()
    
    def plot_regime_distributions(self, save_path: str = None) -> None:
        """Plot return distributions for each regime"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Plot distributions
        for regime, ax in zip([1, 2], [ax1, ax2]):
            regime_returns = self.returns[self.regimes == regime]
            sns.histplot(regime_returns, kde=True, ax=ax)
            ax.set_title(f'Regime {regime} Returns Distribution')
            ax.set_xlabel('Returns')
            ax.set_ylabel('Frequency')
            
            # Add statistics to plot
            stats_text = f'Mean: {np.mean(regime_returns):.2f}\n'
            stats_text += f'Std: {np.std(regime_returns):.2f}\n'
            stats_text += f'Skew: {stats.skew(regime_returns):.2f}\n'
            stats_text += f'Kurt: {stats.kurtosis(regime_returns):.2f}'
            ax.text(0.95, 0.95, stats_text,
                   transform=ax.transAxes,
                   verticalalignment='top',
                   horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        plt.close()

    def plot_qq_plots(self, save_path: str = None) -> None:
        """Create Q-Q plots for each regime to test for normality"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        for regime, ax in zip([1, 2], [ax1, ax2]):
            regime_returns = self.returns[self.regimes == regime]
            stats.probplot(regime_returns, dist="norm", plot=ax)
            ax.set_title(f'Regime {regime} Q-Q Plot')
        
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        plt.close()

    def perform_statistical_tests(self) -> Dict:
        """Perform detailed statistical tests on the regimes"""
        tests = {}
        
        # Perform tests for each regime
        for regime in [1, 2]:
            regime_returns = self.returns[self.regimes == regime]
            
            # Normality tests
            shapiro_test = stats.shapiro(regime_returns)
            ks_test = stats.kstest(regime_returns, 'norm')
            
            # Stationarity test (Augmented Dickey-Fuller)
            from statsmodels.tsa.stattools import adfuller
            adf_test = adfuller(regime_returns)
            
            # Autocorrelation test (Ljung-Box)
            from statsmodels.stats.diagnostic import acorr_ljungbox
            lb_test = acorr_ljungbox(regime_returns, lags=[10], return_df=True)
            
            tests[f'regime_{regime}'] = {
                'shapiro_test': {
                    'statistic': shapiro_test[0],
                    'p_value': shapiro_test[1]
                },
                'ks_test': {
                    'statistic': ks_test[0],
                    'p_value': ks_test[1]
                },
                'adf_test': {
                    'statistic': adf_test[0],
                    'p_value': adf_test[1]
                },
                'ljung_box_test': {
                    'statistic': lb_test['lb_stat'].values[0],
                    'p_value': lb_test['lb_pvalue'].values[0]
                }
            }
        
        # Test for regime differences
        regime1_returns = self.returns[self.regimes == 1]
        regime2_returns = self.returns[self.regimes == 2]
        
        # Mann-Whitney U test for distribution differences
        mw_test = stats.mannwhitneyu(regime1_returns, regime2_returns)
        tests['regime_comparison'] = {
            'mann_whitney_u_test': {
                'statistic': mw_test[0],
                'p_value': mw_test[1]
            }
        }
        
        return tests

    def generate_summary_report(self) -> Dict:
        """Generate comprehensive backtesting report"""
        risk_measures = self.calculate_regime_risk_measures()
        violations, expected = self.calculate_violation_ratios()
        kupiec_results = self.kupiec_test(violations)
        
        report = {
            'risk_measures': risk_measures,
            'violation_ratio': np.mean(violations),
            'expected_violations': np.mean(expected),
            'kupiec_test': kupiec_results,
            'total_observations': len(self.returns),
            'regime_proportions': {
                'regime_1': np.mean(self.regimes == 1),  # Base regime
                'regime_2': np.mean(self.regimes == 2)   # Spike regime
            }
        }
        
        return report

    def generate_extended_summary_report(self) -> Dict:
        """Generate a comprehensive summary report including statistical tests"""
        basic_report = self.generate_summary_report()
        statistical_tests = self.perform_statistical_tests()
        
        return {
            **basic_report,
            'statistical_tests': statistical_tests
        }

# Example usage
if __name__ == "__main__":
    # Initialize backtester
    backtester = MarkovRegimeBacktest(confidence_level=0.95)
    
    # Load MATLAB results (adjust path as needed)
    # backtester.load_matlab_results('path_to_matlab_results.mat')
    
    # Generate and print summary report
    # report = backtester.generate_summary_report()
    # print("Backtesting Results:")
    # print(pd.DataFrame(report).to_string())