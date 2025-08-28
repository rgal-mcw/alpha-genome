'''
track_data.TrackData objects store model predictions. They have the following properties (using tdata as an example of a track_data.TrackData object)

* tdata.values : store track predictions as a numpy.ndarray.
* tdata.metadata : stores track metadata as a pandas.DataFrame. For each track in the predicted values, there will be a corresponding row in the track metadata describing its origin.
* tdata.uns : contains additional unstructured metadata as a dict.
'''

'''
From Scratch:
You can create your own track_data.TrackData object from strach by specifying ithe values and metadata manually. THe metadata must contain at least the columns name (the names of the tracks) and strand (the stands of DNA that the tracks are on):
'''
from alphagenome.data import genome
from alphagenome.models import dna_client
import numpy as np
import pandas as pd
from alphagenome.data import track_data

# Array has shape (4,3) -> sequence is length 4 and there are 3 tracks:
values = np.array([[0,1,2], [3,4,5], [6,7,8], [9,10,11]]).astype(
        np.float32
        )

# We have both the positive and negative strand values for track1, while track2 contains unstranded data.
metadata = pd.DataFrame({
    'name': ['track1', 'track1', 'track2'],
    'strand': ['+', '-', '.'],
    })
tdata = track_data.TrackData(values = values, metadata=metadata)

'''
Resolution
It's always useful to specify the resolution of the tracks and the genomic interval that they come from, if you have this inforamtion available:
'''

interval = genome.Interval(chromosome='chr1', start=1_000, end=1_004)

tdata = track_data.TrackData(
        values=values, metadata=metadata, resolution=1, interval=interval
        )

# Note that the length of the values has to match up with the interval width and resolution. Here is an example specifying that the values actaully represent 128bp resolution tracks (i.e., each number is a summary over 128 base pairs of DNA):

interval = genome.Interval(chromosome='chr1', start=1_000, end=1_512)

tdata = track_data.TrackData(
        values = values, metadata = metadata, resolution=128, interval=interval
        )

'''
Converting between Resolutions
We can also interconvert between resolution. For example, given 1bp resolution predictions, we can downsample the resolution (by summing adjacent values) and return a sequence of length 2:
'''

interval = genome.Interval(chromosome='chr1', start=1_000, end=1_004)

tdata = track_data.TrackData(
        values = values, metadata=metadata, resolution = 1, interval=interval
        )

tdata = tdata.change_resolution(resolution=2)
print(tdata.values)

# We can also upsample track data to get back to 1bp resolution and a sequence of length 4 by repeating values while preserving the sum:

tdata = tdata.change_resolution(resolution=1)
print(tdata.values)


print("\n --- Filtering --- \n")
'''
Filtering
track_data.TrackData objects can be filtered by the type of DNA strand the tracks are on:
'''
print(
        'Positive strand tracks:',
        tdata.filter_to_positive_strand().metadata.name.values,
        )
print(
        'Negative strand tracks:',
        tdata.filter_to_negative_strand().metadata.name.values,
        )
print('Unstranded tracks:', tdata.filter_to_unstranded().metadata.name.values)

print("\n --- Resizing --- \n")
'''
Resizing
We can resize the `track_data.TrackData` to be either smaller (by cropping):
'''
# Re-instantiating the original trackdata.
tdata = track_data.TrackData(
        values = values, metadata = metadata, resolution = 1, interval = interval
        )

# Resize from width (sequence length) of 4 down to 2.
print(tdata.resize(width=2).values)

# Or bigger (by padding with zeroes)
print(tdata.resize(width=8).values)

print("\n --- Slicing --- \n")
'''
Slicing
We can slice into specific positions of the track_data.TrackData:
'''
# Get the final 2 positions only.
print('slice by position:', tdata.slice_by_positions(start=2, end=4).values)
# Same, but using slice_interval:
print(
        'slice by interval',
        tdata.slice_by_interval(
            genome.Interval(chromosome='chr1', start=1_002, end=1_004)
            ).values,
        )

print("\n --- Subsetting Tracks --- \n")
'''
Subsetting Tracks
subset (and reorder) to specific track names:
'''
# Get only tracks with the name 'track1'.
track1_tdata = tdata.select_tracks_by_name(names='track1')
print(track1_tdata.values)

# The metadata gets automatically filtered to track1 too:
print(track1_tdata.metadata.name.values)

# Finally, if we pass in a stranded genome.Interval or leave unspecified as None when constructing a track_data.TrackDatam, we can reverse complement transform our track values in a strand-aware manner:

interval = genome.Interval(
        chromosome='chr1',
        start=1_000, 
        end=1_004,
        strand='+'
        )

tdata = track_data.TrackData(
        values=values, metadata=metadata, resolution=1, interval=interval
        )

print(tdata.reverse_complement().values)
