Here we provide the ability to upload your own datasets for analysis. The format of these datasets must follow the same pattern as the two datasets we built the app off of.

### Variable Names

You must have the following variable names, which are case sensitive:
* __ID__: This identifies each subject for which data is collected. Subjects could include patients, nights, gene sequences, etc.
* __state__: This identifies the state the subject is in during a time period. States could include types of care, sleep stages, DNA bases, etc.
* __duration__: This is the amount of time the subject spent in a state. This should be integer values. The units are arbitrary, as long as they are integers. Durations could represent time spent in care, time spent in a sleep state, number of identical DNA bases in a row, etc.
* __start_time__: This is the time that the subject entered the state. This should be in the same units as duration.

### Data structure:

Data should be uploaded as a __.csv__ file with the correct variable names. Each row in this data file should represent a period of time that a subject spent in a state. For example, if your row is 


<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:black;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:black;}
.tg .tg-c3ow{border-color:inherit;text-align:center;vertical-align:top}
</style>
<table class="tg" align = "center">
  <tr>
    <th class="tg-c3ow">ID</th>
    <th class="tg-c3ow">state</th>
    <th class="tg-c3ow">duration</th>
    <th class="tg-c3ow">start_time</th>
  </tr>
  <tr>
    <td class="tg-c3ow">1</td>
    <td class="tg-c3ow">"A"</td>
    <td class="tg-c3ow">5</td>
    <td class="tg-c3ow">0</td>
  </tr>
</table>

<br/>

Then subject 1 spent five units of time starting at time 0 in state A. A larger example of a valid data set is given below:
