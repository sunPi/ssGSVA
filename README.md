# ssGSVA
This is just a simple wrapper for the single sample method of the GSVA package. It offers a docopt header, so it can be run via the console. The packages should be automatically installed, if necessary please do install any other dependencies. <br>

---

Usage: ssGSEA_quickstart.R [options] <br>

<pre> Options: <br>
          -h --help                          Show this screen. <br>
          --matrix=&lt;string&gt;                  Specify the path to your gene expression matrix (counts, FPKM, etc). <br>
          --gene_signature=&lt;string&gt;          Path to the gene signature(s). <br>
          --verbose=&lt;value&gt;                  If set to T prints all messages [default: F]. <br>
          --version                          ssGSVA Wrapper 1.0. </pre><be>


Copyright [2024] [Jan Rogel]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
