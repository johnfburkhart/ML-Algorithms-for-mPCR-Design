MeasurePoissonLoss = R6::R6Class("MeasurePoissonLoss",
                             inherit = mlr3::MeasureRegr,
                             public = list(
                               initialize = function() {
                                 super$initialize(
                                   # custom id for the measure
                                   id = "poisson_loss",
                                   
                                   # additional packages required to calculate this measure
                                   packages = character(),
                                   
                                   # properties, see below
                                   properties = character(),
                                   
                                   # required predict type of the learner
                                   predict_type = "response",
                                   
                                   # feasible range of values
                                   range = c(-Inf, Inf),
                                   
                                   # minimize during tuning?
                                   minimize = TRUE
                                 )
                               }
                             ),
                             
                             private = list(
                               # custom scoring function operating on the prediction object
                               .score = function(prediction, ...) {
                                 poisson_loss = function(truth, response) {
                                   mean(response - truth*log(response))
                                   
                                 }
                                 
                                 poisson_loss(prediction$truth, prediction$response)
                               }
                             )
)

mlr3::mlr_measures$add("poisson_loss", MeasurePoissonLoss)
