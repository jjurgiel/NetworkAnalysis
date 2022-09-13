# NetworkAnalysis
This set of functions performed Network Contingency Analysis (NCA) on source-level EEG data, and is translated to work with EEG data from the fMRI implemented introduced in _Sripada, C., Kessler, D., Fang, Y., Welsh, R. C., Prem Kumar, K., & Angstadt, M. (2014). Disrupted network architecture of the resting brain in attention‐deficit/hyperactivity disorder. Human brain mapping, 35(9), 4693-4705._

The main goal of NCA is do determine whether a particular brain network, or pair of brain networks, are communicating aberrantly in a particular diagnostic group.

In brief:
* Effective (causal) connectivity is first calculated between pairs of dipoles using the groupSIFT toolbox (https://github.com/sccn/groupSIFT)
* ICA-localized dipoles are projected to a set of 286 parcels described in:
  * _Gordon, E. M., Laumann, T. O., Adeyemo, B., Huckins, J. F., Kelley, W. M., & Petersen, S. E. (2016). Generation and evaluation of a cortical area parcellation from resting-state correlations. Cerebral cortex, 26(1), 288-303._
* These parcels are segregated into 12 classical brain networks (as described in Gordon et al.), e.g. Default Mode, FrontoParietal, Visual, Ventral Attention, etc.
  * For example, the parcel mapping for the Default Mode network is shown below:
  <img src="https://user-images.githubusercontent.com/12466792/189801666-f855f12b-9025-438e-9fcf-d5cdcace5ee0.png" width="300" height="250">
* Network-level statistics (using either t-test, ANOVA, regression model, etc.) are estimated on groupings of between-network and within-network connections, using a permutation-based approach to control for multiple comparisons within networks, in order to determine which networks show above-chance numbers of atypical connections
* False-discovery rate correction is additionally implemented to control for multiple comparisons across networks

Example outcome:

<img src="https://user-images.githubusercontent.com/12466792/189802401-9f382e65-5564-44b2-976d-e64bbcf25f60.png" width="400" height="400">

Here, a 2-group comparison analysis in network connectivity was performed. Red and blue dots indicate the direction of significant differences between groups for a particular connection (e.g., blue dots show that Group 1 > Group 2 for a connection, while red dots show Group 2 > Group 1).

Network cells (i.e., bigger squares) are shaded based on survival of permutation-based control, where color shading represents the percent of significant connections within that cell which are significant in a particular direction (i.e., if all dots are blue, meaning Group 1 > Group 2, then the cell will be shaded blue because Group 1 > Group 2 at the network level).



### NOTICE:
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
