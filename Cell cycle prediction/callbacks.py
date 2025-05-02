import os
import torch
import numpy as np
import cv2
import pickle
from tensorboardX import SummaryWriter
from utils import bytescale
import torch.nn.functional as F
from copy import deepcopy
from imageio import imwrite
from matplotlib import pyplot as plt

class Callback(object):
    """
    Abstract base class used to build new callbacks.
    """

    def __init__(self):
        self.trainer = None
        self.estimator = None
        self.metrics_collection = None
        self.current_results = None
        self.training = True

    def set_trainer(self, trainer):
        self.trainer = trainer
        self.metrics_collection = trainer.metrics_collection
        self.estimator = trainer.estimator
        self.current_results = trainer.current_results

    def on_batch_begin(self, batch):
        pass

    def on_batch_end(self, batch):
        pass

    def on_epoch_begin(self, epoch):
        pass

    def on_epoch_end(self, epoch):
        pass

    def on_train_begin(self):
        pass

    def on_train_end(self):
        pass


class Callbacks(Callback):
    def __init__(self, callbacks):
        super().__init__()
        if isinstance(callbacks, Callbacks):
            callbacks = callbacks.callbacks
        self.callbacks = callbacks
        if callbacks is None:
            self.callbacks = []

    def set_trainer(self, trainer):
        for callback in self.callbacks:
            callback.set_trainer(trainer)

    def on_batch_begin(self, batch, training):
        for callback in self.callbacks:
            callback.training = training
            callback.on_batch_begin(batch)

    def on_batch_end(self, batch):
        for callback in self.callbacks:
            callback.on_batch_end(batch)

    def on_epoch_begin(self, epoch):
        for callback in self.callbacks:
            callback.on_epoch_begin(epoch)

    def on_epoch_end(self, epoch):
        for callback in self.callbacks:
            callback.on_epoch_end(epoch)

    def on_train_begin(self):
        for callback in self.callbacks:
            callback.on_train_begin()

    def on_train_end(self):
        for callback in self.callbacks:
            callback.on_train_end()


def extract_model(model):
    while True:
        if isinstance(model, torch.nn.DataParallel):
            model = model.module
        else:
            return model


class ModelSaver(Callback):
    def __init__(self, save_every, save_name, best_only=True):
        super().__init__()
        self.save_every = save_every
        self.save_name = save_name
        self.best_only = best_only

    def on_epoch_end(self, epoch):
        loss = float(self.metrics_collection.val_metrics["loss"])
        need_save = not self.best_only
        if epoch % self.save_every == 0:
            if loss < self.metrics_collection.best_loss:
                self.metrics_collection.best_loss = loss
                self.metrics_collection.best_epoch = epoch
                need_save = True

            if need_save:
                torch.save(
                    deepcopy(extract_model(self.estimator.model)),
                    os.path.join(self.estimator.save_path, self.save_name).format(
                        epoch=epoch, loss="{:.2}".format(loss)
                    ),
                )


def save_checkpoint(epoch, model_state_dict, optimizer_state_dict, path):
    torch.save(
        {
            "epoch": epoch + 1,
            "state_dict": model_state_dict,
            "optimizer": optimizer_state_dict,
        },
        path,
    )


class CheckpointSaver(Callback):
    def __init__(self, save_every, save_name):
        super().__init__()
        self.save_every = save_every
        self.save_name = save_name

    def on_epoch_end(self, epoch):
        loss = float(self.metrics_collection.val_metrics["loss"])
        if epoch % self.save_every == 0:
            save_checkpoint(
                epoch,
                extract_model(self.estimator.model).state_dict(),
                self.estimator.optimizer.state_dict(),
                os.path.join(self.estimator.save_path, self.save_name).format(
                    epoch=epoch, loss="{:.2}".format(loss)
                ),
            )


class ModelFreezer(Callback):
    def on_epoch_begin(self, epoch):
        warmup = self.estimator.config.warmup
        if hasattr(self.estimator.model, "encoder_stages"):
            encoder_stages = self.estimator.model.encoder_stages
        elif hasattr(self.estimator.model, "module") and hasattr(
            self.estimator.model.module, "encoder_stages"
        ):
            # if the model is an instance of nn.DataParallel
            encoder_stages = self.estimator.model.module.encoder_stages
        else:
            print("WARNING: model freezer is not working.")
            return
        if epoch < warmup:
            for p in encoder_stages.parameters():
                p.requires_grad = False
            if self.estimator.config.input_channel_num != 3:
                for p in encoder_stages[0][0].parameters():
                    p.requires_grad = True

            # for param_group in self.estimator.optimizer.param_groups:
            #     param_group['lr'] = 1e-5
        if epoch == warmup:
            for p in encoder_stages.parameters():
                p.requires_grad = True


class TensorBoard(Callback):
    def __init__(self, logdir):
        super().__init__()
        self.logdir = logdir
        self.writer = None

    def on_train_begin(self):
        os.makedirs(self.logdir, exist_ok=True)
        self.writer = SummaryWriter(self.logdir)

    def on_epoch_end(self, epoch):
        for k, v in self.metrics_collection.train_metrics.items():
            self.writer.add_scalar("train/{}".format(k), float(v), global_step=epoch)

        for k, v in self.metrics_collection.val_metrics.items():
            self.writer.add_scalar("val/{}".format(k), float(v), global_step=epoch)

        for idx, param_group in enumerate(self.estimator.optimizer.param_groups):
            lr = param_group["lr"]
            self.writer.add_scalar(
                "group{}/lr".format(idx), float(lr), global_step=epoch
            )

    def on_train_end(self):
        self.writer.close()


class PredictionSaver(Callback):
    def __init__(self, output_path, scale_target=False, report_interval=10):
        super().__init__()
        self.report_interval = report_interval
        self.epoch = 0
        self.training_batches = []
        self.validation_batches = []
        self.output_path = output_path
        self.scale_target = scale_target

    def on_epoch_begin(self, epoch):
        self.epoch = epoch
        self.training_batches = []
        self.validation_batches = []

    def on_epoch_end(self, epoch):
        print(self.metrics_collection.val_metrics)
        if not self.training:
            self.save_results(self.validation_batches, prefix="Validation")
            self.save_results(self.training_batches, prefix="Training")
        else:
            self.save_results(self.training_batches, prefix="Training")

    def save_results(self, batches, prefix):
        for bk, batch_results in enumerate(batches):
            #print("keys: ", batch_results.keys())
            for k in batch_results.keys():
                batch = batch_results[k]
                if type(batch) is list:
                    continue
                data = batch.cpu().numpy()
                for b in range(data.shape[0]): 
                    if k == "inputs":
                        im = bytescale(data[b]).astype("uint8")
                    else:
                        if self.scale_target:
                            im = data[b]
                            imx = im.max()
                            if imx != 0:
                                im = im / imx
                        else:
                            im = data[b]
                        if im.ndim == 3:
                            im = np.clip(im * 255, 0, 255).astype("uint8")
                        else:
                            im = np.clip(im * 255, 0, 255).astype("uint8")

                    if im.ndim == 3:
                        im = im.transpose(1, 2, 0)
            
                    if k == "outputs":
                        print(im.shape)
                        print('Output prediction range:', im.shape, im[:,:,0].max(),im[:,:,1].max(),im[:,:,2].min(), im[:,:,2].max())
                    if self.epoch % self.report_interval == 0:
                        os.makedirs(
                            os.path.join(
                                self.output_path, "epoch_{}".format(self.epoch)
                            ),
                            exist_ok=True,
                        )
                        imwrite(
                            os.path.join(
                                self.output_path,
                                "epoch_{}/{}_{}_epoch_{}_#batch_{}_#item_{}_{}.png".format(
                                    self.epoch,
                                    prefix,
                                    batch_results.get("image_name")[b],
                                    self.epoch,
                                    bk,
                                    b,
                                    k,
                                ),
                            ),
                            im,
                        )
                    else:
                        imwrite(
                            os.path.join(
                                self.output_path,
                                "{}_#batch_{}_#item_{}_{}.png".format(prefix, bk, b, k),
                            ),
                            im,
                        )

    def on_batch_end(self, batch):
        if self.epoch % self.report_interval == 0 and not self.training:
            self.validation_batches.append(
                {k: v for k, v in self.current_results.items()}
            )
        else:
            self.training_batches = [{k: v for k, v in self.current_results.items()}]


class HistorySaver(Callback):
    def __init__(self, logs_dir, output_path, report_interval=10):
        super().__init__()
        self.logs_dir = logs_dir
        self.output_path = output_path
        if report_interval < 10:
            self.report_interval = 10
        else:
            self.report_interval = report_interval
        self.epochs = []
        self.history_training = []
        self.history_validation = []
        self.best_epoch = None
        self.best_loss = None

    def on_epoch_begin(self, epoch):
        if os.path.exists(os.path.join(self.logs_dir, "history.pkl")):
            with open(os.path.join(self.logs_dir, "history.pkl"),"rb") as f:
                history_dict = pickle.load(f)
                f.close()
            self.epochs = history_dict['epochs']
            self.history_training = history_dict['train_loss']
            self.history_validation = history_dict['val_loss']
            self.best_loss = min(history_dict['val_loss'])
            self.best_epoch = history_dict['epochs'][history_dict['val_loss'].index(self.best_loss)]
        
        self.epochs.append(epoch)

    def on_epoch_end(self, epoch):
        self.history_training.append(
            float(self.metrics_collection.train_metrics["loss"])
        )
        self.history_validation.append(
            float(self.metrics_collection.val_metrics["loss"])
        )
        if self.best_loss == None or self.metrics_collection.best_loss < self.best_loss:
            self.best_loss = self.metrics_collection.best_loss
            self.best_epoch = self.metrics_collection.best_epoch

        history_dict = {'epochs': self.epochs, 'train_loss': self.history_training, 'val_loss': self.history_validation}
        with open(os.path.join(self.logs_dir, "history.pkl"),"wb") as f:
            pickle.dump(history_dict, f)
            f.close()

        if epoch % self.report_interval == 0:
            self.save_history(
                self.epochs[-self.report_interval :],
                self.history_training[-self.report_interval :],
                self.history_validation[-self.report_interval :],
                title="Short_history.png",
            )

        self.save_history(
            self.epochs,
            self.history_training,
            self.history_validation,
            best_epoch=self.best_epoch,
            best_loss=self.best_loss,
            title="History.png",
        )

    def save_history(
        self,
        epochs,
        train_losses,
        val_losses,
        best_epoch=None,
        best_loss=None,
        title=None,
    ):
        plt.figure()
        plt.title(title)
        plt.plot(epochs, train_losses, color="blue", label="Training loss")
        plt.plot(epochs, val_losses, color="darkorange", label="Validation loss")
        if (best_epoch != None) and (best_loss != None):
            plt.scatter(best_epoch, best_loss, s=100)
            plt.text(best_epoch, best_loss, "best_loss = " + "{:.1e}".format(best_loss))
        if max(epochs) > self.report_interval:
            plt.xticks(
                np.arange(
                    min(epochs),
                    max(epochs) + 1,
                    (max(epochs) + 1 - min(epochs)) // self.report_interval,
                )
            )
        else:
            plt.xticks(epochs)
        plt.yticks(
            np.arange(
                0,
                np.max(val_losses + train_losses),
                np.max(val_losses + train_losses) / 10,
            )
        )
        plt.legend(loc="best", shadow=True)

        plt.xlabel("Epoch")
        plt.ylabel("Loss")
        plt.savefig(os.path.join(self.output_path, title), bbox_inches="tight")
        plt.close()
