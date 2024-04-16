```python
from datetime import datetime

date = datetime.now().strftime("%Y%m%d at %I.%M.%S %p")
print(f"Last Saved: {date}\n\n")

!jupyter nbconvert --to markdown $PWD/02_python_data_vizualization.ipynb --output-dir $PWD/notebooks
!cp $PWD/notebooks/02_python_data_vizualization.md $PWD/02_python_data_visualization.Rmd

```

    Last Saved: 20240415 at 01.27.01 PM
    
    
    [NbConvertApp] Converting notebook /Users/tmbmacbookair/Dropbox/2024/TMB2024/COMPLETED/imap-data-exploration/02_python_data_vizualization.ipynb to markdown
    [NbConvertApp] Writing 446 bytes to /Users/tmbmacbookair/Dropbox/2024/TMB2024/COMPLETED/imap-data-exploration/notebooks/02_python_data_vizualization.md



```python

```
