import React, { useEffect, useMemo, useRef, useState } from 'react';
import {
  Box,
  Button,
  createStyles,
  FormControl,
  FormControlLabel,
  FormGroup,
  Input,
  LinearProgress,
  Radio,
  RadioGroup,
  Switch,
  TextField,
  Theme,
  Typography,
  withStyles,
  WithStyles,
} from '@material-ui/core';
import { grey } from '@material-ui/core/colors';
import { useHistory } from 'react-router-dom';

import { faCaretDown, faChartBar } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';

import { useDispatch, useSelector } from 'react-redux';
import bbox from '@turf/bbox';

import DatePicker from 'react-datepicker';

import { extractPropsFromURL } from './util';
import {
  getBoundaryLayerSingleton,
  LayerDefinitions,
} from '../../../config/utils';
import {
  AggregationOperations,
  BoundaryLayerProps,
  NSOLayerProps,
  WMSLayerProps,
  isStatistic,
  LayerKey,
} from '../../../config/types';
import { LayerData } from '../../../context/layers/layer-data';
import { layerDataSelector } from '../../../context/mapStateSlice/selectors';
import { Extent } from '../Layers/raster-utils';
import { availableDatesSelector } from '../../../context/serverStateSlice';
import {
  AnalysisDispatchParams,
  analysisResultSelector,
  clearAnalysisResult,
  isAnalysisLayerActiveSelector,
  isAnalysisLoadingSelector,
  requestAndStoreAnalysis,
  setIsMapLayerActive,
} from '../../../context/analysisResultStateSlice';
import AnalysisTable from './AnalysisTable';
import {
  getAnalysisTableColumns,
  downloadCSVFromTableData,
} from '../../../utils/analysis-utils';
import LayerDropdown from './LayerDropdown';
import { ActionTypes, useAnalyserReducer } from './AnalyserReducer';

const boundaryLayer = getBoundaryLayerSingleton();

function Analyser({ classes }: AnalyserProps) {
  const dispatch = useDispatch();
  const history = useHistory();
  const boundaryLayerData = useSelector(layerDataSelector(boundaryLayer.id)) as
    | LayerData<BoundaryLayerProps>
    | undefined;

  const availableDates = useSelector(availableDatesSelector);
  const analysisResult = useSelector(analysisResultSelector);

  const isAnalysisLoading = useSelector(isAnalysisLoadingSelector);
  const isMapLayerActive = useSelector(isAnalysisLayerActiveSelector);

  const [isAnalyserFormOpen, setIsAnalyserFormOpen] = useState(false);
  const [isTableViewOpen, setIsTableViewOpen] = useState(true);
  const activeUrl: React.MutableRefObject<string> = useRef<string>('');

  // form elements
  const [form, setForm] = useAnalyserReducer();

  // set default date after dates finish loading and when hazard layer changes
  useEffect(() => {
    const dates =
      form.hazardLayerId !== 'placeholder'
        ? availableDates[
            (LayerDefinitions[form.hazardLayerId] as WMSLayerProps)
              ?.serverLayerName
          ]
        : null;
    if (!dates || dates.length === 0) {
      setForm({ type: ActionTypes.SET_SELECTED_DATE, date: null });
    } else {
      setForm({
        type: ActionTypes.SET_SELECTED_DATE,
        date: dates[dates.length - 1],
      });
    }
  }, [availableDates, form.hazardLayerId, setForm]);

  const adminBoundariesExtent = useMemo(() => {
    if (!boundaryLayerData) {
      // not loaded yet. Should be loaded in MapView
      return null;
    }
    return bbox(boundaryLayerData.data) as Extent; // we get extents of admin boundaries to give to the api.
  }, [boundaryLayerData]);

  const statisticOptions = Object.entries(AggregationOperations).map(stat => (
    <FormControlLabel
      key={stat[0]}
      value={stat[1]}
      control={
        <Radio className={classes.radioOptions} color="default" size="small" />
      }
      label={stat[0]}
    />
  ));

  const clearAnalysis = () => {
    dispatch(clearAnalysisResult());

    setForm({ type: ActionTypes.CLEAR_FORM });

    // Reset URL in browser.
    history.replace('');
  };

  const runAnalyser = async () => {
    if (!adminBoundariesExtent) {
      return;
    } // hasn't been calculated yet

    // Detect from issues and throw errors.
    if (!form.selectedDate) {
      throw new Error('Date must be given to run analysis');
    }

    if (
      form.hazardLayerId === 'placeholder' ||
      form.baselineLayerId === 'placeholder'
    ) {
      throw new Error('Layer should be selected to run analysis');
    }

    // Build URL string.
    const analyserShareURL: string = Object.keys(form)
      .map((prop: string) => {
        return [prop, form[prop]].map(encodeURIComponent).join('=');
      })
      .join('&');

    // Make sure to run teh nalyser onoly when needed.
    if (activeUrl.current !== analyserShareURL) {
      // Store the new URL.
      activeUrl.current = analyserShareURL;

      const selectedHazardLayer = LayerDefinitions[
        form.hazardLayerId
      ] as WMSLayerProps;
      const selectedBaselineLayer = LayerDefinitions[
        form.baselineLayerId
      ] as NSOLayerProps;

      const params: AnalysisDispatchParams = {
        hazardLayer: selectedHazardLayer,
        baselineLayer: selectedBaselineLayer,
        date: form.selectedDate,
        statistic: form.statistic,
        extent: adminBoundariesExtent,
        threshold: {
          above: parseFloat(form.aboveThreshold) || undefined,
          below: parseFloat(form.belowThreshold) || undefined,
        },
      };

      await dispatch(requestAndStoreAnalysis(params));

      // Set share URL into the browser.
      history.replace(`?share=true&${analyserShareURL}`);
    }
  };

  useEffect(() => {
    // Early return if data not loaded.
    if (!boundaryLayerData || !availableDates) {
      return;
    }

    const {
      hazardLayerParamId,
      baselineLayerParamId,
      selectedParamDate,
      statisticParam,
      aboveThresholdParam,
      belowThresholdParam,
      fromURL,
    } = extractPropsFromURL(history.location.search);

    // Set data from the URL.
    if (fromURL) {
      // Correct the statistic if invalid.
      let statisticCorrected: AggregationOperations =
        AggregationOperations.Mean;

      if (isStatistic(statisticParam)) {
        // eslint-disable-next-line fp/no-mutation
        statisticCorrected = statisticParam as AggregationOperations;
      }

      // Set all from data from the URL.
      setForm({
        type: ActionTypes.SET_FORM,
        params: {
          hazardLayerId: hazardLayerParamId as LayerKey,
          baselineLayerId: baselineLayerParamId as LayerKey,
          statistic: statisticCorrected,
          selectedDate: selectedParamDate,
          belowThreshold: belowThresholdParam,
          aboveThreshold: aboveThresholdParam,
        },
      });

      // Avoid Running Analyser if required data are invalid or it's not a share link.
      if (
        form.hazardLayerId !== 'placeholder' &&
        form.baselineLayerId !== 'placeholder' &&
        form.selectedDate &&
        fromURL
      ) {
        runAnalyser();
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [
    boundaryLayerData,
    availableDates,
    history.location.search,
    form.selectedDate,
  ]);

  return (
    <div className={classes.analyser}>
      <Button
        variant="contained"
        color="primary"
        onClick={() => {
          setIsAnalyserFormOpen(!isAnalyserFormOpen);
        }}
      >
        <FontAwesomeIcon
          style={{ marginRight: '10px', fontSize: '1.6em' }}
          icon={faChartBar}
        />
        <Typography variant="body2">Run Analysis</Typography>
        <FontAwesomeIcon icon={faCaretDown} style={{ marginLeft: '10px' }} />
      </Button>

      <Box
        className={classes.analyserMenu}
        width={isAnalyserFormOpen ? 'min-content' : 0}
        padding={isAnalyserFormOpen ? '10px' : 0}
      >
        {isAnalyserFormOpen ? (
          <div>
            <div className={classes.newAnalyserContainer}>
              <div className={classes.analyserOptions}>
                <Typography variant="body2">Hazard Layer</Typography>
                <LayerDropdown
                  type="wms"
                  value={form.hazardLayerId}
                  setValue={(paramLayerKey: LayerKey) => {
                    setForm({
                      type: ActionTypes.SET_HAZARD_LAYER_ID,
                      layerKey: paramLayerKey,
                    });
                  }}
                  title="Hazard Layer"
                  className={classes.selector}
                  placeholder="Choose hazard layer"
                />
              </div>
              <div className={classes.analyserOptions}>
                <Typography variant="body2">Statistic</Typography>
                <FormControl component="div">
                  <RadioGroup
                    name="statistics"
                    value={form.statistic}
                    onChange={(
                      paramEvent: React.ChangeEvent<HTMLInputElement>,
                    ) => {
                      setForm({
                        type: ActionTypes.SET_STATISTIC,
                        statistic: paramEvent.target
                          .value as AggregationOperations,
                      });
                    }}
                    row
                  >
                    {statisticOptions}
                  </RadioGroup>
                </FormControl>
              </div>
              <div className={classes.analyserOptions}>
                <Typography variant="body2">Baseline Layer</Typography>
                <LayerDropdown
                  type="nso"
                  value={form.baselineLayerId}
                  setValue={(paramLayerKey: LayerKey) => {
                    setForm({
                      type: ActionTypes.SET_BASELINE_LAYER_ID,
                      layerKey: paramLayerKey,
                    });
                  }}
                  title="Baseline Layer"
                  className={classes.selector}
                  placeholder="Choose baseline layer"
                />
              </div>
              <div className={classes.analyserOptions}>
                <Typography variant="body2">Threshold</Typography>
                <TextField
                  id="filled-number"
                  error={!!form.thresholdError}
                  helperText={form.thresholdError}
                  className={classes.numberField}
                  label="Min"
                  type="number"
                  value={form.aboveThreshold}
                  onChange={(
                    paramEvent: React.ChangeEvent<HTMLInputElement>,
                  ) => {
                    setForm({
                      type: ActionTypes.SET_ABOVE_THRESHOLD,
                      value: paramEvent.target.value,
                    });
                  }}
                  variant="filled"
                />
                <TextField
                  id="filled-number"
                  label="Max"
                  className={classes.numberField}
                  value={form.belowThreshold}
                  onChange={(
                    paramEvent: React.ChangeEvent<HTMLInputElement>,
                  ) => {
                    setForm({
                      type: ActionTypes.SET_BELOW_THRESHOLD,
                      value: paramEvent.target.value,
                    });
                  }}
                  type="number"
                  variant="filled"
                />
              </div>
              <div className={classes.analyserOptions}>
                <Typography variant="body2">Date</Typography>
                <DatePicker
                  selected={
                    form.selectedDate ? new Date(form.selectedDate) : null
                  }
                  onChange={date => {
                    setForm({
                      type: ActionTypes.SET_SELECTED_DATE,
                      date: date?.getTime() || form.selectedDate,
                    });
                  }}
                  maxDate={new Date()}
                  todayButton="Today"
                  peekNextMonth
                  showMonthDropdown
                  showYearDropdown
                  dropdownMode="select"
                  customInput={<Input />}
                  popperClassName={classes.calendarPopper}
                  includeDates={
                    form.hazardLayerId !== 'placeholder'
                      ? availableDates[
                          (LayerDefinitions[
                            form.hazardLayerId
                          ] as WMSLayerProps).serverLayerName
                        ]?.map(d => new Date(d)) || []
                      : []
                  }
                />
              </div>
            </div>

            {!isAnalysisLoading && analysisResult && (
              <>
                <FormGroup>
                  <FormControlLabel
                    control={
                      <Switch
                        color="default"
                        checked={isMapLayerActive}
                        onChange={e =>
                          dispatch(setIsMapLayerActive(e.target.checked))
                        }
                      />
                    }
                    label="Map View"
                  />
                  <FormControlLabel
                    control={
                      <Switch
                        color="default"
                        checked={isTableViewOpen}
                        onChange={e => setIsTableViewOpen(e.target.checked)}
                      />
                    }
                    label="Table View"
                  />
                </FormGroup>
                {isTableViewOpen && (
                  <AnalysisTable
                    tableData={analysisResult.tableData}
                    columns={getAnalysisTableColumns(analysisResult)}
                  />
                )}
                <Button
                  className={classes.innerAnalysisButton}
                  onClick={() => downloadCSVFromTableData(analysisResult)}
                >
                  <Typography variant="body2">Download</Typography>
                </Button>
                <Button
                  className={classes.innerAnalysisButton}
                  onClick={clearAnalysis}
                >
                  <Typography variant="body2">Clear Analysis</Typography>
                </Button>
              </>
            )}
            {!analysisResult && (
              <Button
                className={classes.innerAnalysisButton}
                onClick={runAnalyser}
                disabled={
                  !!form.thresholdError || // if there is a threshold error
                  !form.selectedDate || // or date hasn't been selected
                  !form.hazardLayerId || // or hazard layer hasn't been selected
                  !form.baselineLayerId || // or baseline layer hasn't been selected
                  isAnalysisLoading // or analysis is currently loading
                }
              >
                <Typography variant="body2">Run Analysis</Typography>
              </Button>
            )}
            {isAnalysisLoading ? <LinearProgress /> : null}
          </div>
        ) : null}
      </Box>
    </div>
  );
}

const styles = (theme: Theme) =>
  createStyles({
    analyser: {
      zIndex: theme.zIndex.drawer,
      position: 'absolute',
      top: 2,
      left: 2,
      textAlign: 'left',
    },
    analyserMenu: {
      backgroundColor: '#5A686C',
      maxWidth: '100vw',
      color: 'white',
      overflowX: 'hidden',
      whiteSpace: 'nowrap',
      borderTopRightRadius: '10px',
      borderBottomRightRadius: '10px',
      height: 'auto',
      maxHeight: '60vh',
    },
    analyserButton: {
      height: '36px',
      'margin-left': '3px',
    },
    analyserOptions: {
      padding: '5px 0px',
    },
    newAnalyserContainer: {
      padding: '5px',
      marginTop: '10px',
    },
    radioOptions: {
      '&.Mui-checked': { color: grey[50] },
      padding: '2px 10px 2px 20px',
    },
    innerAnalysisButton: {
      backgroundColor: '#3d474a',
      margin: '10px',
      '&.Mui-disabled': { opacity: 0.5 },
    },
    selectorLabel: {
      '&.Mui-focused': { color: 'white' },
    },
    selector: {
      margin: '5px',
    },
    numberField: {
      paddingLeft: '10px',
      marginTop: '10px',
      width: '85.5px',
      '& .Mui-focused': { color: 'white' },
    },
    calendarPopper: {
      zIndex: 3,
    },
  });

interface AnalyserProps extends WithStyles<typeof styles> {}

export default withStyles(styles)(Analyser);
